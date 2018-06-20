import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog
import os
import sys


Labels=['Vol','Shell','AClust','Pick','TotVol']

class DataRecordForm(tk.Frame):
    """The input form for our widgets"""

    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.dirname=os.path.expanduser('~')

        # A dict to keep track of input widgets
        self.inputs = {}
        self.inputs['Input']=tk.StringVar()
        self.inputs['Output']=tk.StringVar()
        self.inputs['Begin']=tk.DoubleVar()
        self.inputs['End']=tk.DoubleVar()
        self.dir=tk.StringVar()
        self.dir.set(self.dirname)
        self.relInput=tk.StringVar()
        self.relOutput=tk.StringVar()

        self.vcmd = (self.register(self.onValidate),'%S')

        fromWhere = tk.LabelFrame(self, text="Trajectory")

        myDir=ttk.Label(fromWhere,text='Working Directory:')
        myDir.grid(row=0,column=1, columnspan=1,pady=(1,4))

        myDirLabel=ttk.Label(fromWhere,textvariable=self.dir)
        myDirLabel.grid(row=1,column=1,columnspan=1,pady=(5,20))

        myInput=ttk.Button(fromWhere,text='Input',command=self.on_Input)
        myInputLabel=ttk.Label(fromWhere,textvariable=self.relInput)
        myInput.grid(row=2,column=0, columnspan=1,pady=(5,10))
        myInputLabel.grid(row=2,column=1, columnspan=1,pady=(5,10))

        myOutput=ttk.Button(fromWhere,text='Output',command=self.on_Output)
        myOutputLabel=ttk.Label(fromWhere,textvariable=self.relOutput)
        myOutput.grid(row=2,column=2, columnspan=1,pady=(5,10))        
        myOutputLabel.grid(row=2,column=3, columnspan=1,pady=(5,10))

        self.makeEntry(fromWhere,text='Begin',row=3,column=0)
        self.makeEntry(fromWhere,text='End',row=3,column=2)        
        fromWhere.grid(row=3,column=0,columnspan=2)

        myOptions = tk.LabelFrame(self, text="Options")
        ttk.Label(myOptions,text='Pick an option').grid(row=0,column=0,pady=(0,10),padx=(5,1))
        self.which=tk.StringVar()
        self.inputs['Which']=ttk.Combobox(myOptions, textvariable=self.which, values=Labels,width=10)
        self.inputs['Which'].grid(row=0,column=1,pady=(0,10),padx=(5,1))
        self.whichArg=tk.StringVar()
        self.makeEntry(myOptions,text='Argument',row=0,column=2)
        myOptions.grid(row=4,column=0,sticky=tk.W)
        style = ttk.Style()
        style.configure("RW.TLabel", foreground="red")
        self.error_message = ttk.Label(parent, text='', style="RW.TLabel")
        self.error_message.grid(row=5, pady=(0,20))

    def getDir(self):
        self.dir.set(self.dirname)
        return self.dir

    def omy(self,event):
        messagebox.showinfo('Baby', 'Message content')

    def makeEntry(self,parent,text=None,textvar=None,row=0,column=0,width=10):
        ttk.Label(parent,text=text).grid(row=row,column=column,pady=(0,10),padx=(5,1))
        self.inputs[text]=ttk.Entry(parent,text=text,width=width,validate="key", validatecommand=self.vcmd)
        self.inputs[text].grid(row=row,column=column+1,pady=(0,10),padx=5)
        
    def get(self):
        data = {}
        for key in self.inputs:
            data[key] = self.inputs[key].get()
        return data

    def reset(self):
        for widget in self.inputs.values():
            widget.set('')
    def on_Dir(self):
        dir=filedialog.askdirectory(initialdir = os.path.expanduser('~'),title = "Select dircetory")
        self.dirname=dir

    def on_Input(self):
        filename =  filedialog.askopenfilename(initialdir = os.path.expanduser('~'),title = "Select file",filetypes=(("json files", "*.json"),
                                           ("All files", "*.*") ))
        self.dirname=os.path.dirname(filename)
        self.dir.set(self.dirname)
        var=os.path.relpath(filename,self.dirname)
        self.relInput.set(var)
        self.inputs['Input'].set(filename)

    def on_Output(self):
        filename =  filedialog.asksaveasfile(initialdir = self.dirname,title = "Save to",filetypes=(("output files", "*.dat"),
                                           ("All files", "*.*") )).name
        var=os.path.relpath(filename,self.dirname)
        self.relOutput.set(var)
        self.inputs['Output'].set(filename)


    def onValidate(self,S):
        self.error_message.config(text='')            
        if S.isdigit() or S == '.':
            return True
        else:
            self.error_message.config(text='Invalid character %s \n name can only have alphabets and spaces' % S)            
            return False
