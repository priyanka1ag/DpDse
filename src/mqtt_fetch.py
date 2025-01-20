import numpy as np
import paho.mqtt.client as mqtt
from io import StringIO
import pandas as pd
import uuid
import time
from datetime import datetime
import threading


# Global variable to hold the received data
m = None

# Callback when the client connects to the broker
def on_connect_sub(client, userdata, flags, rc):
    if rc == 0:
        print("SUB: Connected successfully!")
    else:
        print(f"SUB: Failed to connect, return code {rc}")



# Callback function for handling incoming MQTT messages
def on_message(client, userdata, msg):
    #print("message is: ", msg.payload.decode())
    global m
    # Convert the message payload from bytes to a string, then convert to a numpy array
    try:
        m = pd.read_json(StringIO(msg.payload.decode()))
        m = pd.DataFrame(m)
        print("Dataframe: ", m)
        # Update measurements (u & z)
        # predict
        # correct
        print("time is: ", time.time()) 
        # upon receiving measurements execute predict and correct steps - to update both u and z
        userdata.predict()
        userdata.correct()

    except ValueError:
        print("Failed to convert message payload to a dataframe.")



topic = "network_10_nodes"
client_id = str(uuid.uuid4())
print(f"Using unique Client ID: { client_id  }")
broker_ip = "broker.hivemq.com"
port = 1883
'''
cl = mqtt.Client(mqtt.CallbackAPIVersion.VERSION1, client_id = client_id)
cl.on_connect = on_connect_sub
cl.on_message = on_message
cl.enable_logger()
cl.connect("broker.hivemq.com", 1883)
time.sleep(2)
cl.loop_start()
#cl.loop_forever()

cl.subscribe(topic)
time.sleep(2)
cl.loop_stop()
'''




def storedata_once(client1, topic):
    while True:
        try:
            client1.subscribe(topic)
        except:
            print("Unexpected error:", + " ts: " + str(datetime.now()))
        else:
            break

def run_predict(stop_event, dse):
    while not stop_event.is_set():
        time.sleep(0.001)
        dse.predict()

def storedata_repeatedly(topic, dse):
    sec = 0
    client_id = str(uuid.uuid4())
    client1 = mqtt.Client(mqtt.CallbackAPIVersion.VERSION1, client_id = client_id, userdata=dse) # for v2.x version of Paho-mqtt
    client1.on_message = on_message
    client1.on_connect = on_connect_sub
    # Start continuous function in a separate thread
    stop_event = threading.Event()
    thread = threading.Thread(target=run_predict, args=(stop_event, dse))
    thread.start()
    client1.connect(broker_ip, port)
    client1.loop_start()
    time.sleep(0.2)
    while sec<60:
        storedata_once(client1, topic)
        time.sleep(0.1)
        sec = sec + 1
    return


def rerun_every_minute(topic, dse):
    while 1:
        storedata_repeatedly(topic, dse)
        print("next run")


#rerun_every_minute(topic)