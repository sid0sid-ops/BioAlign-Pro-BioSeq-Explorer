const puppeteer = require('puppeteer-core');
const fs = require('fs');
const path = require('path');
const http = require('http');

const CHROME_PATH = 'C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe';
const APP_URL = 'http://127.0.0.1:3838';

function wait(ms) {
    return new Promise(resolve => setTimeout(resolve, ms));
}

// Function to poll the server until it is ready
function waitForServer(url, timeoutMs = 60000) {
    return new Promise((resolve, reject) => {
        const startTime = Date.now();
        const interval = setInterval(() => {
            if (Date.now() - startTime > timeoutMs) {
                clearInterval(interval);
                reject(new Error(`Server at ${url} did not respond within ${timeoutMs}ms`));
                return;
            }
            
            http.get(url, (res) => {
                if (res.statusCode === 200) {
                    clearInterval(interval);
                    resolve();
                }
            }).on('error', () => {
                // Ignore error and retry
            });
        }, 1000);
    });
}

async function run() {
    console.log(`Polling for Shiny app at ${APP_URL}...`);
    try {
        await waitForServer(APP_URL);
        console.log('Shiny app is online!');
    } catch (err) {
        console.error(err.message);
        process.exit(1);
    }

    console.log('Launching Chrome...');
    const browser = await puppeteer.launch({
        executablePath: CHROME_PATH,
        headless: true,
        args: ['--no-sandbox', '--disable-setuid-sandbox']
    });

    try {
        const page = await browser.newPage();
        await page.setViewport({ width: 1400, height: 900 });

        console.log(`Navigating to ${APP_URL}...`);
        await page.goto(APP_URL, { waitUntil: 'networkidle2' });

        // Ensure directory Images exists
        const imgDir = path.join(__dirname, '..', 'Images');
        if (!fs.existsSync(imgDir)) {
            fs.mkdirSync(imgDir, { recursive: true });
        }

        console.log('Waiting for dashboard cards...');
        await page.waitForSelector('#tab_manager-qa_card_codon_usage', { timeout: 10000 });
        await wait(2000); // Wait for animations and initial rendering

        // --- 1. CODON USAGE TAB ---
        console.log('Opening Codon Usage tab...');
        await page.click('#tab_manager-qa_card_codon_usage');
        await page.waitForSelector('#tab_manager-codon_usage-codon_usage_tool_root', { timeout: 15000 });
        await wait(1500);

        console.log('Running Codon Analysis...');
        await page.waitForSelector('#tab_manager-codon_usage-run_full', { timeout: 5000 });
        await page.click('#tab_manager-codon_usage-run_full');
        
        // Wait for analysis to finish (wait 8 seconds)
        console.log('Waiting for Codon analysis to finish...');
        await wait(8000);

        console.log('Taking codon_overview.png screenshot...');
        await page.screenshot({ path: path.join(imgDir, 'codon_overview.png') });

        console.log('Navigating to Codon Bias tab...');
        const biasTabSelector = '#tab_manager-codon_usage-result_tabs a[data-value="Codon Bias"]';
        await page.waitForSelector(biasTabSelector, { timeout: 5000 });
        await page.click(biasTabSelector);
        await wait(3000); // wait for Plotly/Echarts to load
        console.log('Taking codon_bias.png screenshot...');
        await page.screenshot({ path: path.join(imgDir, 'codon_bias.png') });

        console.log('Navigating to Optimization Studio...');
        const optTabSelector = '#tab_manager-codon_usage-result_tabs a[data-value="Optimization Studio"]';
        await page.waitForSelector(optTabSelector, { timeout: 5000 });
        await page.click(optTabSelector);
        await wait(1500);

        console.log('Running Codon Optimization...');
        const optBtnSelector = '#tab_manager-codon_usage-run_optimization';
        await page.waitForSelector(optBtnSelector, { timeout: 5000 });
        await page.click(optBtnSelector);
        await wait(6000); // Wait for optimization calculation and diff table
        console.log('Taking codon_optimization.png screenshot...');
        await page.screenshot({ path: path.join(imgDir, 'codon_optimization.png') });

        // --- 2. MOTIF SEARCH TAB ---
        console.log('Navigating back to Dashboard...');
        const dashTabSelector = '#tab_manager-workspace_tabs a[data-value="Dashboard"]';
        await page.waitForSelector(dashTabSelector, { timeout: 5000 });
        await page.click(dashTabSelector);
        await wait(2000);

        console.log('Opening Motif Search tab...');
        await page.click('#tab_manager-qa_card_motif_search');
        await page.waitForSelector('#tab_manager-motif_search-motif_root', { timeout: 15000 });
        await wait(1500);

        // Click Run Analysis in Motif Search
        console.log('Running Motif Search...');
        const motifSearchBtn = '#tab_manager-motif_search-btn_search';
        await page.waitForSelector(motifSearchBtn, { timeout: 5000 });
        await page.click(motifSearchBtn);
        await wait(8000); // Wait for search and charts to load

        console.log('Taking motif_overview.png screenshot...');
        await page.screenshot({ path: path.join(imgDir, 'motif_overview.png') });

        console.log('Navigating to Motif Visualizations...');
        const motifVisTab = '#tab_manager-motif_search-result_tabs a[data-value="visualizations"]';
        await page.waitForSelector(motifVisTab, { timeout: 5000 });
        await page.click(motifVisTab);
        await wait(4000);
        console.log('Taking motif_visualizations.png screenshot...');
        await page.screenshot({ path: path.join(imgDir, 'motif_visualizations.png') });

        console.log('Navigating to Motif Structure-Aware Analysis...');
        const motifStrTab = '#tab_manager-motif_search-result_tabs a[data-value="structure"]';
        await page.waitForSelector(motifStrTab, { timeout: 5000 });
        await page.click(motifStrTab);
        await wait(4000);
        console.log('Taking motif_structure.png screenshot...');
        await page.screenshot({ path: path.join(imgDir, 'motif_structure.png') });

        console.log('All screenshots captured successfully!');
    } catch (err) {
        console.error('Automation failed:', err);
    } finally {
        await browser.close();
        console.log('Browser closed.');
    }
}

run();
