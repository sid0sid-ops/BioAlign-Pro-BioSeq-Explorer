packages = ['easyocr', 'cv2', 'numpy', 'matplotlib', 'pandas', 'openpyxl', 'pytesseract', 'tesseract']
for p in packages:
    try:
        __import__(p)
        print(f"Package {p} is INSTALLED")
    except ImportError:
        print(f"Package {p} is NOT installed")
