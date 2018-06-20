import Voro
from subprocess import *
import tkinter as tk
from tkinter import ttk
from time import *



root=tk.Tk()
canvas=tk.Text(root,width=80, height=24)
canvas.grid(row=0,column=0)
proc=Popen('/Users/marchi/git/openTRJ/python/Voro.py', stdout=PIPE,stderr=PIPE)
proc=proc.communicate()
canvas.insert(tk.END,proc)
sleep(1)
proc=Popen('/Users/marchi/git/openTRJ/python/Voro.py', stdout=PIPE,stderr=PIPE)
proc=proc.communicate()
canvas.insert(tk.END,proc)
sleep(1)
proc=Popen('/Users/marchi/git/openTRJ/python/Voro.py', stdout=PIPE,stderr=PIPE)
proc=proc.communicate()
canvas.insert(tk.END,proc)

root.mainloop()