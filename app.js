// TCGA-PRAD Research Explorer Application Logic

document.addEventListener('DOMContentLoaded', async () => {
  // Navigation Tabs
  const navButtons = document.querySelectorAll('.nav-item');
  const panels = document.querySelectorAll('.view-panel');

  navButtons.forEach(btn => {
    btn.addEventListener('click', () => {
      const targetId = btn.getAttribute('data-tab');
      
      navButtons.forEach(b => b.classList.remove('active'));
      panels.forEach(p => p.classList.remove('active'));

      btn.classList.add('active');
      const targetPanel = document.getElementById(`tab-${targetId}`);
      if (targetPanel) {
        targetPanel.classList.add('active');
      }
    });
  });

  // Load Verified Research Artifact
  try {
    const res = await fetch('static/research/tcga_prad_metrics.json');
    if (res.ok) {
      const data = await res.json();
      console.log('Loaded verified research metrics:', data);
    }
  } catch (err) {
    console.warn('Artifact static load notice:', err);
  }
});
