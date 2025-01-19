import matplotlib
import numpy as np
import customtkinter as ctk
import splines as spl
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

matplotlib.use('TkAgg')
ctk.set_appearance_mode("light")
ctk.set_default_color_theme("green")


class BaseFrame(ctk.CTkFrame):
    def __init__(self, master, f, fder, f2der, a, b, **kwargs):
        super().__init__(master, **kwargs)
        self.f = f
        self.fder = fder
        self.f2der = f2der
        self.a = a
        self.b = b

        # надписи
        self.label = ctk.CTkLabel(self, text="Введите число отрезков разбиения")
        self.label.place(x=35, y=15)

        self.n_label = ctk.CTkLabel(self, text="n = ")
        self.n_label.place(x=35, y=45)
        self.n_var = ctk.IntVar(self, value=256)

        # Поля
        self.n_entry = ctk.CTkEntry(self, width=75, textvariable=self.n_var)
        self.n_entry.place(x=65, y=45)

        # Кнопки
        self.run_button = ctk.CTkButton(self, text="Рассчитать", command=self.run)
        self.run_button.place(x=35, y=85)

        # Таблица
        self.place_tables()

        # Графики
        self.figure = Figure(figsize=(7, 7.5), dpi=100, facecolor="#cfcfcf")
        self.figure_canvas = FigureCanvasTkAgg(self.figure, self)
        self.axes = self.figure.subplots(nrows=2, ncols=1, sharex=True)
        self.figure_canvas.get_tk_widget().place(x=830, y=75)

        # Текстовые поля
        self.textbox = ctk.CTkTextbox(self, width=250, height=100, border_color="#2cc985", border_width=4)

        # Выбор графиков
        self.segm_gr = ctk.CTkSegmentedButton(master=self, values=["Функции", "Производные",
                                                                   "2 Производные", 'Погрешности'],
                                              command=self.segm_gr)
        self.segm_gr.place(x=730, y=65)

    def segm_gr(self, value):
        if value == "Функции":
            self.plot_fun()
        elif value == "Производные":
            self.plot_der()
        elif value == "2 Производные":
            self.plot_2der()
        else:
            self.plot_errors()

    def run(self):
        self.splines = spl.find_splines(self.n_var.get(), self.f, self.a, self.b, 0, 0)
        self.N = 2 * self.n_var.get()
        self.X = np.linspace(self.a, self.b, self.N+1)
        self.F = np.array([self.f(xi) for xi in self.X])
        self.Fder = np.array([self.fder(xi) for xi in self.X])
        self.F2der = np.array([self.f2der(xi) for xi in self.X])
        self.S = []
        self.Sder = []
        self.S2der = []
        i = 0
        for spline in self.splines:
            while self.X[i] < spline.x[1]:
                fun = spline.get_func()
                funder = spline.get_func_derivative()
                fun2der = spline.get_func_2derivative()
                self.S.append(fun(self.X[i]))
                self.Sder.append(funder(self.X[i]))
                self.S2der.append(fun2der(self.X[i]))
                i += 1
        self.S.append(self.F[-1])
        self.Sder.append(self.Fder[-1])
        self.S2der.append(self.F2der[-1])
        self.S = np.array(self.S)
        self.Sder = np.array(self.Sder)
        self.S2der = np.array(self.S2der)

        self.error = np.absolute(self.F - self.S)
        self.der_error = np.absolute(self.Fder - self.Sder)
        self.sec_der_error = np.absolute(self.F2der - self.S2der)

        self.show_table()
        self.show_info()

    def plot_fun(self):
        self.axes[0].cla()
        self.axes[0].set_title('График сплайна:')
        self.axes[0].set_ylabel("S(x)")
        self.axes[0].set_xlabel("x")
        for spline in self.splines:
            x = np.linspace(spline.x[0], spline.x[1], 5)
            f = spline.get_func()
            y = [f(xi) for xi in x]
            self.axes[0].plot(x, y, c='g', linewidth=2)
        self.axes[0].relim()

        self.axes[1].cla()
        self.axes[1].set_title('График функции:')
        self.axes[1].set_ylabel("f(x)")
        self.axes[1].set_xlabel("x")
        x = np.linspace(self.a, self.b, 101)
        y = np.array([self.f(xi) for xi in x])
        self.axes[1].plot(x, y, c='r', linewidth=2)
        self.axes[1].relim()

        self.figure_canvas.draw()

    def plot_der(self):
        self.axes[0].cla()
        self.axes[0].set_title('График производной сплайна:')
        self.axes[0].set_ylabel("S'(x)")
        self.axes[0].set_xlabel("x")
        for spline in self.splines:
            x = np.linspace(spline.x[0], spline.x[1], 5)
            f = spline.get_func_derivative()
            y = [f(xi) for xi in x]
            self.axes[0].plot(x, y, c='g', linewidth=2)
        self.axes[0].relim()

        self.axes[1].cla()
        self.axes[1].set_title('График производной функции:')
        self.axes[1].set_ylabel("f'(x)")
        self.axes[1].set_xlabel("x")
        x = np.linspace(self.a, self.b, 51)
        y = np.array([self.fder(xi) for xi in x])
        self.axes[1].plot(x, y, c='r', linewidth=2)
        self.axes[1].relim()

        self.figure_canvas.draw()

    def plot_2der(self):
        self.axes[0].cla()
        self.axes[0].set_title('График 2 производной сплайна:')
        self.axes[0].set_ylabel("S''(x)")
        self.axes[0].set_xlabel("x")
        for spline in self.splines:
            x = np.linspace(spline.x[0], spline.x[1], 5)
            f = spline.get_func_2derivative()
            y = [f(xi) for xi in x]
            self.axes[0].plot(x, y, c='g', linewidth=2)
        self.axes[0].relim()

        self.axes[1].cla()
        self.axes[1].set_title('График 2 производной функции:')
        self.axes[1].set_ylabel("f''(x)")
        self.axes[1].set_xlabel("x")
        x = np.linspace(self.a, self.b, 51)
        y = np.array([self.f2der(xi) for xi in x])
        self.axes[1].plot(x, y, c='r', linewidth=2)
        self.axes[1].relim()

        self.figure_canvas.draw()

    def plot_errors(self):
        self.axes[0].cla()
        self.axes[0].set_title('График погрешности:')
        self.axes[0].set_ylabel("F(x) - S(x)")
        self.axes[0].set_xlabel("x")
        self.axes[0].plot(self.X, self.error, c='r', linewidth=2)
        self.axes[0].relim()

        self.axes[1].cla()
        self.axes[1].set_title('График погрешности производных:')
        self.axes[1].set_xlabel("x")
        self.axes[1].plot(self.X, self.der_error, c='g', linewidth=2, label='Погр. произв.')
        self.axes[1].plot(self.X, self.sec_der_error, c='b', linewidth=2, label='Погр. 2 произв.')
        self.axes[1].legend()
        self.axes[1].relim()

        self.figure_canvas.draw()

    def show_table(self):
        self.table.destroy()
        self.scrollbar_y.destroy()
        self.table_comp.destroy()
        self.comp_scrollbar_y.destroy()

        self.place_tables()

        # таблица сплайнов
        self.table["column"] = ['i', 'x(i-1)', 'xi', 'ai', 'bi', 'ci', 'di']
        self.table["show"] = "headings"
        for column in self.table["columns"]:
            self.table.heading(column, text=column)
        self.table.column(0, width=30)
        for i, spline in enumerate(self.splines):
            self.table.insert("", "end", values=(i + 1, spline.x[0],
                              spline.x[1], spline.a, spline.b, spline.c, spline.d))

        # таблица сравнения
        self.table_comp["column"] = ['j', 'xj', 'F(xj)', 'S(xj)', 'F(xj) - S(xj)',
                                     "F'(xj)", "S'(xj)", "F'(xj) - S'(xj)",
                                     "F''(xj)", "S''(xj)", "F''(xj) - S''(xj)"]
        self.table_comp["show"] = "headings"
        for column in self.table_comp["columns"]:
            self.table_comp.heading(column, text=column)
        self.table_comp.column(0, width=30)
        for j in range(self.N + 1):
            self.table_comp.insert("", "end",
                                   values=(j, self.X[j], self.F[j], self.S[j], self.F[j] - self.S[j],
                                           self.Fder[j], self.Sder[j], self.Fder[j] - self.Sder[j],
                                           self.F2der[j], self.S2der[j], self.F2der[j] - self.S2der[j]))

    def place_tables(self):
        self.table = ttk.Treeview(self, height=20)
        self.scrollbar_y = ttk.Scrollbar(self, orient=ctk.VERTICAL, command=self.table.yview)
        self.table.place(x=35, y=165, width=740, height=265)
        self.table.configure(yscroll=self.scrollbar_y.set)
        self.scrollbar_y.place(x=775, y=165, height=265)

        self.table_comp = ttk.Treeview(self, height=20)
        self.comp_scrollbar_y = ttk.Scrollbar(self, orient=ctk.VERTICAL, command=self.table_comp.yview)
        self.table_comp.place(x=35, y=480, width=740, height=265)
        self.table_comp.configure(yscroll=self.comp_scrollbar_y.set)
        self.comp_scrollbar_y.place(x=775, y=480, height=265)

    def show_info(self):
        self.textbox.destroy()
        self.textbox = ctk.CTkTextbox(self, width=350, height=100, border_color="#2cc985", border_width=4)
        self.textbox.insert("0.0", f"Сетка сплайна n = {self.n_var.get()}\n"
                            f"Контрольная сетка N = {self.N}\n"
                            f"Погрешность сплайна на контрольной сетке\n"
                            f"max|F(xj) - S(xj)| = {np.max(self.error)},\n"
                            f"при x = {self.X[np.argmax(self.error)]}\n"
                            f"Погрешность производной на контрольной сетке\n"
                            f"max|F'(xj) - S'(xj)| = {np.max(self.der_error)},\n"
                            f"при x = {self.X[np.argmax(self.der_error)]}\n"
                            f"Погрешность второй производной на контрольной сетке\n"
                            f"max|F''(xj) - S''(xj)| = {np.max(self.sec_der_error)},\n"
                            f"при x = {self.X[np.argmax(self.sec_der_error)]}\n")
        self.textbox.place(x=290, y=15)


class TestTaskFrame(ctk.CTkFrame):
    def __init__(self, master, **kwargs):
        super().__init__(master, **kwargs)

        self.frame = BaseFrame(self, spl.phi, spl.phi_derivative, spl.phi_2derivative, -1, 1, width=1250, height=670)
        self.frame.place(x=0, y=0)


class MainTask1Frame(ctk.CTkFrame):
    def __init__(self, master, **kwargs):
        super().__init__(master, **kwargs)
        # надписи
        self.fun_label = ctk.CTkLabel(self, text='Выберите функцию: ')
        self.fun_label.place(x=15, y=5)

        # выбор функции
        self.fun_segm = ctk.CTkSegmentedButton(self, values=['ln(x+1)/(x+1)', 'ln(x+1)/x', 'sin(x+1)/x'],
                                               command=self.fun_selection)
        self.fun_segm.place(x=150, y=5)

        # фрейм
        self.frame = BaseFrame(self, spl.f1, spl.f1_der, spl.f1_2der, 0.2, 2, width=1250, height=620)
        self.frame.place(x=0, y=45)

    def fun_selection(self, value):
        if value == 'ln(x+1)/(x+1)':
            self.f = spl.f1
            self.fder = spl.f1_der
            self.f2der = spl.f1_2der
            self.a = 0.2
            self.b = 2
        elif value == 'ln(x+1)/x':
            self.f = spl.f2
            self.fder = spl.f2_der
            self.f2der = spl.f2_2der
            self.a = 2
            self.b = 4
        else:
            self.f = spl.f3
            self.fder = spl.f3_der
            self.f2der = spl.f3_2der
            self.a = 1
            self.b = np.pi

        self.frame.place_forget()
        self.frame = BaseFrame(self, self.f, self.fder, self.f2der, self.a, self.b, width=1250, height=630)
        self.frame.place(x=0, y=45)


class MainTask2Frame(ctk.CTkFrame):
    def __init__(self, master, **kwargs):
        super().__init__(master, **kwargs)
        # надписи
        self.fun_label = ctk.CTkLabel(self, text='Выберите функцию: ')
        self.fun_label.place(x=15, y=5)

        self.cos_label = ctk.CTkLabel(self, text='+ cos(10x)')
        self.cos_label.place(x=385, y=5)

        # выбор функции
        self.fun_segm = ctk.CTkSegmentedButton(self, values=['ln(x+1)/(x+1)', 'ln(x+1)/x', 'sin(x+1)/x'],
                                               command=self.fun_selection)
        self.fun_segm.place(x=150, y=5)

        # фрейм
        self.frame = BaseFrame(self, spl.f1_osc, spl.f1_osc_der, spl.f1_osc_2der, 0.2, 2, width=1250, height=620)
        self.frame.place(x=0, y=45)

    def fun_selection(self, value):
        if value == 'ln(x+1)/(x+1)':
            self.f = spl.f1_osc
            self.fder = spl.f1_osc_der
            self.f2der = spl.f1_osc_2der
            self.a = 0.2
            self.b = 2
        elif value == 'ln(x+1)/x':
            self.f = spl.f2_osc
            self.fder = spl.f2_osc_der
            self.f2der = spl.f2_osc_2der
            self.a = 2
            self.b = 4
        else:
            self.f = spl.f3_osc
            self.fder = spl.f3_osc_der
            self.f2der = spl.f3_osc_2der
            self.a = 1
            self.b = np.pi

        self.frame.place_forget()
        self.frame = BaseFrame(self, self.f, self.fder, self.f2der, self.a, self.b, width=1250, height=630)
        self.frame.place(x=0, y=45)


class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Интерполяция кубическими сплайнами")
        self.geometry("1280x720")

        # Фреймы
        self.test_task_frame = TestTaskFrame(master=self, width=1250, height=680)
        self.main_task1_frame = MainTask1Frame(master=self, width=1250, height=680)
        self.main_task2_frame = MainTask2Frame(master=self, width=1250, height=680)

        # Надписи
        self.label = ctk.CTkLabel(self, text="Выберите задачу: ")
        self.label.place(x=15, y=5)

        # Кнопки
        self.segm_button = ctk.CTkSegmentedButton(master=self, values=["Тестовая задача",
                                                                       "Основная задача 1",
                                                                       "Основная задача 2"],
                                                  command=self.segmented_button_callback)
        self.segm_button.place(x=140, y=5)

    def segmented_button_callback(self, value):
        if value == "Тестовая задача":
            self.main_task2_frame.place_forget()
            self.main_task1_frame.place_forget()
            self.test_task_frame.place(x=15, y=40)
        elif value == "Основная задача 1":
            self.main_task2_frame.place_forget()
            self.test_task_frame.place_forget()
            self.main_task1_frame.place(x=15, y=40)
        else:
            self.main_task1_frame.place_forget()
            self.test_task_frame.place_forget()
            self.main_task2_frame.place(x=15, y=40)
