import websocket
import time
import threading

def on_message(ws, message):
    print("Received from server:", message)

def on_error(ws, error):
    print("WebSocket Error:", error)

def on_close(ws, close_status_code, close_msg):
    print("WebSocket Closed")

def on_open(ws):
    print("WebSocket Opened")
    # Shiny protocol handshake: Shiny expects client to send its initial input state or at least keep connection alive.
    # Let's wait a bit to let session run and then close.
    def run():
        time.sleep(5)
        print("Closing WebSocket connection...")
        ws.close()
    threading.Thread(target=run).start()

if __name__ == "__main__":
    ws = websocket.WebSocketApp("ws://127.0.0.1:3838/",
                              on_open=on_open,
                              on_message=on_message,
                              on_error=on_error,
                              on_close=on_close)
    ws.run_forever()
