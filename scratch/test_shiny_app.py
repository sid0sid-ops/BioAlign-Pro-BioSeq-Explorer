import urllib.request

try:
    response = urllib.request.urlopen("http://127.0.0.1:3838/")
    print("Status code:", response.getcode())
    print("Headers:", response.info())
    print("Body snippet:", response.read()[:500])
except Exception as e:
    print("Connection failed:", e)
