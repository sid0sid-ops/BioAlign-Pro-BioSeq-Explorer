(function () {
  let timer = null;

  function resizeAllCharts() {
    // Resize all ECharts instances
    if (window.echarts) {
      document.querySelectorAll("[_echarts_instance_]").forEach(function (node) {
        try {
          const instance = echarts.getInstanceByDom(node);
          if (instance) {
            instance.resize();
          }
        } catch (error) {
          console.warn("ECharts resize failed", error);
        }
      });
    }
    // Resize all Plotly instances
    if (window.Plotly) {
      document.querySelectorAll(".js-plotly-plot, .shiny-plotly-output, .plotly").forEach(function (node) {
        try {
          if (node.layout) {
            Plotly.Plots.resize(node);
          }
        } catch (error) {
          console.warn("Plotly resize failed", error);
        }
      });
    }
  }

  function resizeSoon() {
    clearTimeout(timer);
    timer = setTimeout(function () {
      resizeAllCharts();
      // Repeated dispatches to catch transition delays
      setTimeout(resizeAllCharts, 150);
      setTimeout(resizeAllCharts, 350);
      setTimeout(resizeAllCharts, 650);
    }, 120);
  }

  // Bind to browser and Bootstrap events
  window.addEventListener("resize", resizeSoon);
  document.addEventListener("shown.bs.tab", resizeSoon);
  document.addEventListener("shown.bs.modal", resizeSoon);
  document.addEventListener("hidden.bs.modal", resizeSoon);
  
  // Bind to Shiny hooks
  document.addEventListener("shiny:value", resizeSoon);
  document.addEventListener("shiny:bound", resizeSoon);
  document.addEventListener("shiny:idle", resizeSoon);

  // Bind click handler for expand modal clicks and sidebar toggle
  document.addEventListener("click", function (event) {
    if (
      event.target.closest(".expand-view-btn") ||
      event.target.closest("[data-action='expand']") ||
      event.target.closest(".modal") ||
      event.target.closest(".sidebar-toggle") ||
      event.target.closest("#btn_sidebar_toggle") ||
      event.target.closest(".plot-expand-button")
    ) {
      resizeSoon();
    }
  });

  // ResizeObserver for direct container layout transitions
  if ("ResizeObserver" in window) {
    const observer = new ResizeObserver(function (entries) {
      entries.forEach(function (entry) {
        if (entry.contentRect.width > 120) {
          resizeSoon();
        }
      });
    });

    function observeContainers() {
      document.querySelectorAll(".codon-plot-card-v2 .echarts4r, .codon-plot-card-v2 .html-widget, .square-chart-wrap").forEach(function (el) {
        observer.observe(el);
      });
    }

    observeContainers();
    setTimeout(observeContainers, 500);
    setTimeout(observeContainers, 1500);
  }

  setTimeout(resizeSoon, 300);
})();
