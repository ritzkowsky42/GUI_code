import tkinter

root = tkinter.Tk()
root.geometry("600x400")
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


tkinter.Label(root, text='First Name', textvariable=Fname).grid(row=0)
tkinter.Label(root, text='Last Name', textvariable=Lname).grid(row=1)
e1 = tkinter.Entry(root)
e2 = tkinter.Entry(root)
e1.grid(row=0, column=1)
e2.grid(row=1, column=1)
print(Fname, Lname)
tkinter.mainloop()

