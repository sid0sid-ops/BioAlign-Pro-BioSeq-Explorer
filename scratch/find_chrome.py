import os

paths = [
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
    r"C:\Users\Siddharth Tripathi\AppData\Local\Google\Chrome\Application\chrome.exe"
]

found = False
for p in paths:
    if os.path.exists(p):
        print(f"FOUND: {p}")
        found = True
        break

if not found:
    print("Chrome not found in standard paths.")
