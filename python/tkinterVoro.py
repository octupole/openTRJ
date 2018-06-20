import tkinter as tk
from tkinter import ttk

from datetime import datetime
import os
import csv
from DataRecordForm import *
from time import sleep
import Voro
from subprocess import *

# Start coding here

commands={'Input':'-f','Output':'-o','Begin':'-b','End':'-e','Which':'--which'}

class Application(tk.Tk):
    """Application root window"""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.title("Voro Entry Application")

        self.resizable(width=True, height=True)
        mylabel=ttk.Label(self,text="Voro Entry Application",
            font=("TkDefaultFont", 16)
        )
        mylabel.grid(row=0,columnspan=4,pady=(10,10))

        self.recordform = DataRecordForm(self)
        self.recordform.grid(row=1, column=0,padx=40,pady=(40,5))


        self.run = ttk.Button(self, text="Run", command=self.on_running)
        self.run.grid(row=3,column=0,sticky=(tk.E),padx=(0,35),pady=(0,20))
        self.termf = ttk.Frame(self, height=100, width=500)
        self.termf.grid(row=4,column=0,sticky=(tk.E),padx=(0,35),pady=(0,20))
        


    def on_running(self):
        mydata=self.recordform.get()
        for key in mydata:
            if key != 'Input' and key != 'Argument':
                if not mydata[key]:
                    mydata[key]=None
        argv=[]
        for key in mydata:
            if key != 'Which' and key != 'Argument':
                a=[commands[key],mydata[key]]
                if mydata[key]:
                    argv+=a
            elif key == 'Which':
                a=[commands[key]+'='+mydata[key]]
                if mydata['Argument']:
                    a=[commands[key]+'='+mydata[key],mydata['Argument']]
                argv+=a
        self.argv=argv
        output=tk.Canvas(self.termf, width=100, height=100, bg="white")
        proc = Popen(Voro.main(self.argv), stdout=PIPE)
        # proc = proc.communicate()
        # output.insert(tk.END, proc)




if __name__ == "__main__":
    app = Application()
    app.mainloop()
    # app.destroy()
    # Voro.main(app.argv)