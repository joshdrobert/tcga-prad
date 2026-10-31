from __future__ import annotations

import http.server
import socketserver
import threading
from functools import partial
from pathlib import Path
try:
    from playwright.sync_api import sync_playwright
except ImportError:
    sync_playwright = None

WEB_DIR = Path(__file__).resolve().parents[1] / "web"


class QuietSimpleHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format, *args):
        pass


def run_ui_tests():
    port = 8899
    handler = partial(QuietSimpleHTTPRequestHandler, directory=str(WEB_DIR))

    # Allow port reuse
    socketserver.TCPServer.allow_reuse_address = True
    httpd = socketserver.TCPServer(("127.0.0.1", port), handler)
    server_thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    server_thread.start()

    url = f"http://127.0.0.1:{port}/index.html"
    print(f"Testing UI at {url}...")

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)

        viewports = [
            ("Desktop", {"width": 1280, "height": 800}),
            ("Tablet", {"width": 768, "height": 1024}),
            ("Mobile", {"width": 375, "height": 667}),
        ]

        for name, vp in viewports:
            console_errors = []
            page = browser.new_page(viewport=vp)

            page.on("console", lambda msg: console_errors.append(msg.text) if msg.type == "error" else None)
            page.on("pageerror", lambda exc: console_errors.append(str(exc)))

            page.goto(url, wait_until="domcontentloaded")
            page.wait_for_selector(".kpi-value")

            # Check title
            title = page.title()
            assert "TCGA-PRAD" in title, f"Unexpected page title on {name}: {title}"

            # Check KPI values
            assert page.locator(".kpi-value").first.is_visible(), f"KPI cards not visible on {name}"

            # Test tab navigation
            tab_buttons = page.locator(".nav-item").all()
            assert len(tab_buttons) == 8, f"Expected 8 navigation tabs, found {len(tab_buttons)} on {name}"

            for btn in tab_buttons:
                tab_id = btn.get_attribute("data-tab")
                btn.click()
                panel = page.locator(f"#tab-{tab_id}")
                assert panel.is_visible(), f"Panel #{tab_id} not visible after click on {name}"
                assert "active" in (panel.get_attribute("class") or ""), f"Panel #{tab_id} missing active class on {name}"

            # Check images loaded
            images = page.locator("img.chart-img").all()
            for img in images:
                src = img.get_attribute("src")
                assert src, f"Chart image missing src on {name}"

            assert len(console_errors) == 0, f"Found console errors on {name}: {console_errors}"
            print(f"✓ {name} passed with 0 console errors (all 8 tabs active)")
            page.close()

        browser.close()

    httpd.shutdown()
    print("All Playwright UI tests completed successfully!")


def test_ui_playwright():
    if sync_playwright is None:
        import pytest
        pytest.skip("Playwright not installed in current environment")
    run_ui_tests()


if __name__ == "__main__":
    if sync_playwright is None:
        print("Playwright not installed in current environment")
        exit(1)
    run_ui_tests()
