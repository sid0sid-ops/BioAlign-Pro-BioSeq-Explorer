// =====================================================================
// BioSeq Explorer — Custom JavaScript
// =====================================================================

document.addEventListener("DOMContentLoaded", function () {
  // Force light theme – no toggle allowed
  document.body.classList.add("light-theme");

  // Remove theme toggle event listener (disabled)
  // If a toggle exists in the UI, it will be ignored.
  const toggle = document.getElementById("theme_toggle");
  if (toggle) {
    toggle.disabled = true;
    toggle.checked = false; // ensure unchecked (light)
  }

  // ── Sidebar Drawer Toggle (Hide / Show) + Overlay ─────────────────────
  document.addEventListener("click", function (e) {
    const sbToggle = e.target.closest("#btn_sidebar_toggle");
    const overlay = e.target.closest("#sidebar_overlay");
    const sidebar = document.querySelector(".sidebar-left");
    const overlayEl = document.getElementById("sidebar_overlay");
    
    // Toggle via hamburger OR click on overlay to close
    if (sbToggle || overlay) {
      if (sidebar) {
        sidebar.classList.toggle("collapsed");
        const isCollapsed = sidebar.classList.contains("collapsed");
        if (overlayEl) {
          if (isCollapsed) {
            overlayEl.classList.remove("active");
          } else {
            overlayEl.classList.add("active");
          }
        }
        
        // Update ARIA expanded state for accessibility
        const toggleBtn = document.getElementById("btn_sidebar_toggle");
        if (toggleBtn) {
          toggleBtn.setAttribute("aria-expanded", isCollapsed ? "false" : "true");
        }
        
        // Trigger window resize to adjust echarts
        window.dispatchEvent(new Event("resize"));
        setTimeout(() => window.dispatchEvent(new Event("resize")), 150);
        setTimeout(() => window.dispatchEvent(new Event("resize")), 300);
      }
    }
  });


  // ── Adaptive Settings Dropdown ──────────────────────────────────────
  document.addEventListener("click", function(e) {
    const btn = e.target.closest("#btn_settings_dropdown_toggle");
    console.log("Document clicked. Target:", e.target, "Settings Button found:", btn);
    if (btn) {
      const dropdown = btn.closest(".settings-dropdown");
      console.log("Settings Button clicked. Dropdown container:", dropdown);
      if (dropdown) {
        dropdown.classList.toggle("show");
        console.log("Toggled 'show' class. Current classes:", dropdown.className);
      }
    } else {
      // Close dropdown if clicking outside
      const activeDropdown = document.querySelector(".settings-dropdown.show");
      if (activeDropdown && !e.target.closest(".settings-dropdown-content")) {
        activeDropdown.classList.remove("show");
        console.log("Closed active dropdown because click was outside");
      }
    }
  });

  // ── Copy to Clipboard ──────────────────────────────────────────────
  document.addEventListener("click", function (e) {
    const btn = e.target.closest(".btn-copy");
    if (!btn) return;

    const targetId = btn.dataset.target;
    const el = document.getElementById(targetId);
    if (!el) return;

    const text = el.innerText.trim();
    navigator.clipboard.writeText(text).then(function () {
      btn.classList.add("copied");
      const original = btn.innerHTML;
      btn.innerHTML = "✓ Copied";
      setTimeout(function () {
        btn.innerHTML = original;
        btn.classList.remove("copied");
      }, 1500);
    });
  });

  // ── Nucleotide box hover tooltip ───────────────────────────────────
  document.addEventListener("mouseover", function (e) {
    const box = e.target.closest(".nuc-box");
    if (!box) return;
    const nuc = box.textContent.trim();
    const names = { A: "Adenine", T: "Thymine", G: "Guanine", C: "Cytosine" };
    box.title = names[nuc] || nuc;
  });

  // ── Trigger resize when tabs change to fix ECharts squishing and Plotly blurriness ────────
  document.addEventListener("shown.bs.tab", function () {
    window.dispatchEvent(new Event("resize"));
    
    // Explicitly resize Plotly instances to avoid hidden-tab scaling blurriness
    if (window.Plotly) {
      document.querySelectorAll('.plotly .js-plotly-plot').forEach(function(el) {
        try { Plotly.Plots.resize(el); } catch(e) {}
      });
    }

    setTimeout(function() {
      window.dispatchEvent(new Event("resize"));
      if (window.Plotly) {
        document.querySelectorAll('.plotly .js-plotly-plot').forEach(function(el) {
          try { Plotly.Plots.resize(el); } catch(e) {}
        });
      }
    }, 150);

    setTimeout(function() {
      window.dispatchEvent(new Event("resize"));
      if (window.Plotly) {
        document.querySelectorAll('.plotly .js-plotly-plot').forEach(function(el) {
          try { Plotly.Plots.resize(el); } catch(e) {}
        });
      }
    }, 300);
  });

  // ── Strict Accordion Navigation Menus ───────────────────────────────
  document.addEventListener("click", function(e) {
    const btn = e.target.closest(".nav-tab-btn");
    if (!btn) return;
    
    const mode = btn.getAttribute("data-target");
    if(!mode) return;

    const isActive = btn.classList.contains("active");
    
    // Deactivate ALL buttons and collapse ALL panels first (strict accordion)
    const accordion = btn.closest(".yt-accordion-sidebar");
    if(accordion) {
      accordion.querySelectorAll(".nav-tab-btn").forEach(b => b.classList.remove("active"));
      accordion.querySelectorAll(".collapse-panel").forEach(p => {
        p.classList.remove("active-panel");
      });
    }
    
    // If it was NOT active, open it (and it alone)
    if (!isActive) {
      btn.classList.add("active");
      const parentItem = btn.closest(".accordion-item-custom");
      if (parentItem) {
        const targetPanel = parentItem.querySelector(".collapse-panel");
        if (targetPanel) {
          targetPanel.classList.add("active-panel");
        }
      }
      
      // Auto-expand sidebar if it's collapsed
      const sidebar = btn.closest(".sidebar-left");
      if (sidebar && sidebar.classList.contains("collapsed")) {
        sidebar.classList.remove("collapsed");
        const overlayEl = document.getElementById("sidebar_overlay");
        if (overlayEl) overlayEl.classList.add("active");
      }
    }
    
    window.dispatchEvent(new Event("resize"));
  });

  // ── Example Dropdown Item Clicks ──────────────────────────────────
  document.addEventListener("click", function(e) {
    const item = e.target.closest(".example-item");
    if (!item) return;
    e.preventDefault();
    const val = item.getAttribute("data-val");
    // Find closest dropdown button to extract namespace
    const dropdownBtn = item.closest(".dropdown").querySelector(".dropdown-toggle");
    if (!dropdownBtn) return;
    const id = dropdownBtn.id;
    const nsPrefix = id.substring(0, id.indexOf("-") + 1);
    Shiny.setInputValue(nsPrefix + "example_selected", val, {priority: "event"});
  });

  // ── Client-side PNG Export (Full Sequence View) ──────────────────────────
  document.addEventListener("click", function(e) {
    const btn = e.target.closest(".btn-export-png");
    if (!btn) return;
    
    const activePane = document.querySelector(".tab-content > .tab-pane.active");
    if (!activePane) {
      alert("No active sequence view to export.");
      return;
    }
    
    const scrollContainer = activePane.querySelector("div[style*='overflow']");
    const targetEl = scrollContainer || activePane;
    
    const origHeight = targetEl.style.maxHeight;
    const origOverflow = targetEl.style.overflowY;
    
    // Expand to full height so html2canvas captures everything
    targetEl.style.maxHeight = "none";
    targetEl.style.overflowY = "visible";
    
    // Attempt to grab sequence name for the file name
    const titleText = document.querySelector("#seq_name_input")?.value || "BioSeq_Export";
    
    btn.innerHTML = '<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span> Exporting...';
    btn.disabled = true;
    
    setTimeout(() => {
      html2canvas(targetEl, {
        scale: 2,
        useCORS: true,
        backgroundColor: "#ffffff"
      }).then(canvas => {
        targetEl.style.maxHeight = origHeight;
        targetEl.style.overflowY = origOverflow;
        
        btn.innerHTML = '<i class="bi bi-camera"></i> Export to PNG';
        btn.disabled = false;
        
        const link = document.createElement("a");
        link.download = `${titleText}_sequence_view.png`;
        link.href = canvas.toDataURL("image/png");
        link.click();
      }).catch(err => {
        console.error("Export failed:", err);
        targetEl.style.maxHeight = origHeight;
        targetEl.style.overflowY = origOverflow;
        btn.innerHTML = '<i class="bi bi-camera"></i> Export to PNG';
        btn.disabled = false;
        alert("Failed to export PNG. See console for details.");
      });
    }, 100); // small delay to let DOM paint the expanded height
  });

  // ── Client-side Print / PDF Export trigger ──────────────────────────
  document.addEventListener("click", function(e) {
    const btn = e.target.closest(".btn-print-view");
    if (!btn) return;
    window.print();
  });

  // ── Responsive Sequence Viewer Wrapping with ResizeObserver ────────
  var sequenceResizeObserver = null;

  function updateSequenceWrapWidth() {
    var viewer = document.querySelector(".motif-sequence-viewer");
    if (!viewer) return;
    
    // Create a temporary span to measure the exact monospace character width
    var testSpan = document.createElement("span");
    testSpan.style.fontFamily = "'JetBrains Mono', monospace";
    testSpan.style.fontSize = window.getComputedStyle(viewer).fontSize || "14px";
    testSpan.style.letterSpacing = "1.5px";
    testSpan.style.visibility = "hidden";
    testSpan.style.position = "absolute";
    testSpan.style.whiteSpace = "nowrap";
    testSpan.textContent = "A".repeat(100);
    document.body.appendChild(testSpan);
    
    var charWidth = testSpan.offsetWidth / 100;
    document.body.removeChild(testSpan);
    
    var containerWidth = viewer.clientWidth;
    if (containerWidth && charWidth) {
      // 55px line number width + 16px right margin + 28px viewer padding
      var availableWidth = containerWidth - 99; 
      var wrapChars = Math.floor(availableWidth / charWidth);
      wrapChars = Math.max(35, Math.min(250, wrapChars));
      
      // Extract namespace dynamically
      var outputEl = document.querySelector("[id$='results_highlighted_sequence']");
      if (outputEl) {
        var id = outputEl.id;
        var nsPrefix = id.substring(0, id.indexOf("results_highlighted_sequence"));
        var inputId = nsPrefix + "motif_sequence_wrap_width";
        Shiny.setInputValue(inputId, wrapChars);
      }
    }
  }

  function setupSequenceObserver() {
    var viewer = document.querySelector(".motif-sequence-viewer");
    if (!viewer) return;

    if (sequenceResizeObserver) {
      sequenceResizeObserver.disconnect();
    }

    sequenceResizeObserver = new ResizeObserver(function(entries) {
      updateSequenceWrapWidth();
    });

    sequenceResizeObserver.observe(viewer);
  }

  // Bind to resize and tab switches
  window.addEventListener("resize", function() {
    setTimeout(updateSequenceWrapWidth, 50);
  });
  document.addEventListener("shown.bs.tab", function() {
    setTimeout(function() {
      updateSequenceWrapWidth();
      setupSequenceObserver();
    }, 50);
  });
  
  // Also poll on Shiny value updates
  $(document).on("shiny:value", function(event) {
    if (event.target && event.target.id && event.target.id.indexOf("results_highlighted_sequence") !== -1) {
      setTimeout(function() {
        updateSequenceWrapWidth();
        setupSequenceObserver();
      }, 100);
    }
  });

});
