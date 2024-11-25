import tkinter
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,  NavigationToolbar2Tk 
import numpy as np

root = tkinter.Tk()
root.geometry("900x700")
menu = tkinter.Menu(root)
root.config(menu=menu)
filemenu = tkinter.Menu(menu)
menu.add_cascade(label='File', menu=filemenu)
filemenu.add_command(label='New')
filemenu.add_command(label='Open...')
filemenu.add_separator()
filemenu.add_command(label='Exit', command=root.quit)
helpmenu = tkinter.Menu(menu)
menu.add_cascade(label='Help', menu=helpmenu)
helpmenu.add_command(label='About')
Fname = tkinter.StringVar()
Lname = tkinter.StringVar()
Freq = tkinter.DoubleVar()
Phase = tkinter.IntVar()

def submit():

    first_name_str = Fname.get()
    last_name_str = Lname.get()
    
    print("The name is : " + first_name_str + " " + last_name_str)
    
    plot(first_name_str, last_name_str)
    Fname.set("")
    Lname.set("")

  
# plot function is created for  
# plotting the graph in  
# tkinter window 
def plot(first_name_str, last_name_str): 
  
    # the figure that will contain the plot 
    fig = Figure(figsize = (5, 5), 
                 dpi = 100) 
  
    # list of squares 
    t = np.arange(0, 5, 0.01)
    omg = Freq.get()
    phase = Phase.get()
    y = np.cos(t*omg + np.pi*phase/180)
  
    # adding the subplot 
    plot1 = fig.add_subplot(111) 
  
    # plotting the graph 
    plot1.plot(t, y) 
    plot1.set_title("Cos requested by " + first_name_str +" " + last_name_str)
  
    # creating the Tkinter canvas 
    # containing the Matplotlib figure 
    canvas = FigureCanvasTkAgg(fig, 
                               master = root)
    canvas.get_tk_widget().grid(row=7,column=2,columnspan=10,rowspan=10)   
    canvas.draw() 

    toolbarFrame = tkinter.Frame(master=root)
    toolbarFrame.grid(row=18,column=2)
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)
  
  
    # placing the toolbar on the Tkinter window 
    #canvas.grid(row=7, column=2)

tkinter.Label(root, text='First Name').grid(row=0)
tkinter.Label(root, text='Last Name').grid(row=1)
tkinter.Label(root, text='Frequency').grid(row=2)
tkinter.Label(root, text='Phase (deg)').grid(row=3)
e1 = tkinter.Entry(root, textvariable=Fname)
e2 = tkinter.Entry(root, textvariable=Lname)
e3 = tkinter.Entry(root, textvariable=Freq)
e4 = tkinter.Entry(root, textvariable=Phase)
e1.grid(row=0, column=1)
e2.grid(row=1, column=1)
e3.grid(row=2, column=1)
e4.grid(row=3, column=1)

sub_btn=tkinter.Button(root,text = 'Submit', command = submit)
sub_btn.grid(row=5,column=1)

tkinter.mainloop()

