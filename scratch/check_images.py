import os
from PIL import Image

image_dir = r"c:\Users\Siddharth Tripathi\OneDrive\Desktop\bioseq_explorer\Images"
images = sorted([f for f in os.listdir(image_dir) if f.endswith(".png")], key=lambda x: int(x.split(".")[0]) if x.split(".")[0].isdigit() else x)

print("Images in directory:")
for img_name in images:
    img_path = os.path.join(image_dir, img_name)
    try:
        with Image.open(img_path) as img:
            width, height = img.size
            size_kb = os.path.getsize(img_path) / 1024
            print(f"File: {img_name} | Dimensions: {width}x{height} | Size: {size_kb:.2f} KB | Format: {img.format}")
    except Exception as e:
        print(f"Error reading {img_name}: {e}")
