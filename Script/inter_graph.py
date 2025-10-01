import tkinter as tk 
from tkinter import ttk
from tkinter.messagebox import showinfo

studies =[]
def add_study():
    global studies
    name_study=Study.get()
    studies.append(name_study)

def OGQ():

    global selected
    bouton.config(state=tk.DISABLED)
    
    Radio_label=tk.Label(Main_win,text="Est ce que cette étude aura besoin d un Outgroup specific")
    r1 = ttk.Radiobutton(Main_win, text='Yes', value="Y", variable=selected,command=mettre_a_jour)
    r2 = ttk.Radiobutton(Main_win, text='No', value="N", variable=selected,command=mettre_a_jour)
    Radio_label.pack()
    r1.pack()
    r2.pack()




def mettre_a_jour():
    global combien
    global num
    if selected.get() == "Y":  # Si "Oui"
        combien_label.pack()
        num.pack()
    else:  # Si "Non"
        combien_label.pack_forget()
        num.pack_forget()

Main_win=tk.Tk()
selected = tk.StringVar()
stud=tk.StringVar()
inte=tk.IntVar()
label = tk.Label(Main_win, text="Please enter the name of your study and press add")
Study = ttk.Entry(Main_win,textvariable=stud)

bouton = tk.Button(Main_win, text="Add", fg="black",command=OGQ)
combien_label=tk.Label(Main_win,text="Combien d'outgroup aller vous avoir besoin?")
num = ttk.Entry(Main_win,textvariable=inte)
label.pack()
Study.pack()
bouton.pack()





Main_win.mainloop()
print(studies)