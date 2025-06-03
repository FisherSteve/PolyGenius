import math
import tkinter as tk
from tkinter import ttk, messagebox 
from fractions import Fraction 
import numpy as np # Für die allgemeine Nullstellensuche und linspace
import random # Für den "Ich fühle mich glücklich"-Button

# Matplotlib-Integration für Tkinter
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

class CubicPolynomial:
    """
    Repräsentiert ein kubisches Polynom P(x) = ax^3 + bx^2 + cx + d.
    Kann auch zur Darstellung von Polynomen niedrigeren Grades verwendet werden (a=0 etc.).
    Koeffizienten werden intern als Fraction-Objekte gespeichert, um Präzision zu gewährleisten.
    """
    def __init__(self, a, b, c, d):
        try:
            self.a = Fraction(a)
            self.b = Fraction(b)
            self.c = Fraction(c)
            self.d = Fraction(d)
        except (TypeError, ValueError) as e:
            raise ValueError(f"Ungültiger Koeffiziententyp oder -wert: {e}")

    def __str__(self):
        if self.a == 0 and self.b == 0 and self.c == 0 and self.d == 0:
            return "0"
        
        processed_terms = []
        is_first_displayed_term = True
        
        coeff_power_pairs = [
            (self.a, "x^3"),
            (self.b, "x^2"),
            (self.c, "x"),
            (self.d, "")
        ]

        for coeff_frac, power_str in coeff_power_pairs:
            if coeff_frac == 0:
                continue 

            sign_str, val_str = "", ""
            
            if is_first_displayed_term:
                if coeff_frac < 0: sign_str = "-"
            else:
                sign_str = " - " if coeff_frac < 0 else " + "
            
            abs_coeff = abs(coeff_frac)
            
            if abs_coeff.denominator == 1: 
                val_str = str(abs_coeff.numerator) if not (abs_coeff.numerator == 1 and power_str != "") else ""
            else: 
                val_str = f"{abs_coeff.numerator}/{abs_coeff.denominator}"

            if is_first_displayed_term and val_str == "" and sign_str == "" and power_str != "" and coeff_frac > 0 :
                 pass 
            elif is_first_displayed_term and val_str == "" and sign_str == "-" and power_str != "":
                 pass 

            processed_terms.append(f"{sign_str}{val_str}{power_str}")
            is_first_displayed_term = False 
            
        final_str = "".join(processed_terms)
        
        if not final_str: return "0" # Sollte nur passieren, wenn alle Koeffizienten 0 waren
        
        return final_str

    def evaluate(self, x):
        # Für numpy arrays oder einzelne Werte, konvertiere zu float für matplotlib
        if isinstance(x, np.ndarray):
            # Vektorisiere die Fraction-Berechnung oder konvertiere Koeffizienten zu float
            a_f, b_f, c_f, d_f = float(self.a), float(self.b), float(self.c), float(self.d)
            return a_f * x**3 + b_f * x**2 + c_f * x + d_f
        else:
            x_frac = Fraction(x) 
            return self.a * x_frac**3 + self.b * x_frac**2 + self.c * x_frac + self.d

    def derivative1_coeffs(self):
        return [3 * self.a, 2 * self.b, self.c]

    def derivative2_coeffs(self):
        return [6 * self.a, 2 * self.b]

    @classmethod
    def from_saddle_point(cls, x0, D0, P_triple_prime_0, force_integer_coeffs=True):
        f_P_triple_prime_0 = Fraction(P_triple_prime_0)
        f_x0 = Fraction(x0)
        f_D0 = Fraction(D0)
        if f_P_triple_prime_0 == 0:
            raise ValueError("P_triple_prime_0 darf nicht Null sein.")
        if force_integer_coeffs:
            if not (f_P_triple_prime_0.denominator == 1 and f_P_triple_prime_0.numerator % 6 == 0):
                 raise ValueError("Für ganzzahlige Koeffizienten muss P_triple_prime_0 eine ganze Zahl und ein Vielfaches von 6 sein.")
            a_gen = f_P_triple_prime_0 / 6 
        else:
            a_gen = f_P_triple_prime_0 / 6 
        a_p, b_p, c_p, d_p = a_gen, -3*f_x0*a_gen, 3*f_x0**2*a_gen, -f_x0**3*a_gen + f_D0
        if force_integer_coeffs:
            return cls(a_p.limit_denominator(1), b_p.limit_denominator(1), c_p.limit_denominator(1), d_p.limit_denominator(1))
        else:
            return cls(a_p, b_p, c_p, d_p)

    @classmethod
    def from_extrema_locations(cls, q1, q2, K_user, C0):
        f_q1, f_q2, f_K_user, f_C0 = Fraction(q1), Fraction(q2), Fraction(K_user), Fraction(C0)
        if f_q1 == f_q2: raise ValueError("q1 und q2 müssen verschieden sein.")
        if f_K_user == 0: raise ValueError("K_user darf nicht Null sein.")
        a_p,b_p,c_p,d_p = 2*f_K_user, -3*f_K_user*(f_q1+f_q2), 6*f_K_user*f_q1*f_q2, f_C0
        return cls(a_p, b_p, c_p, d_p)

    @classmethod
    def from_extrema_locations_and_values(cls, q1, q2, y1, y2, force_integer_coeffs=False):
        f_q1, f_q2 = Fraction(q1), Fraction(q2)
        f_y1, f_y2 = Fraction(y1), Fraction(y2)
        if f_q1 == f_q2:
            raise ValueError("q1 und q2 (x-Koordinaten der Extrema) müssen verschieden sein.")
        def term_val(x_val, f_q1_in, f_q2_in):
            x = Fraction(x_val)
            return 2*x**3 - 3*(f_q1_in+f_q2_in)*x**2 + 6*f_q1_in*f_q2_in*x
        T1, T2 = term_val(f_q1, f_q1, f_q2), term_val(f_q2, f_q1, f_q2)
        if T1 == T2:
            if f_y1 == f_y2: raise ValueError("K_user unbestimmt: T(q1)=T(q2) und y1=y2.")
            else: raise ValueError("Keine Lösung für K_user: T(q1)=T(q2) aber y1!=y2.")
        K_user = (f_y2 - f_y1) / (T2 - T1)
        C0 = f_y1 - K_user * T1
        
        a_p = 2 * K_user
        b_p = -3 * K_user * (f_q1 + f_q2)
        c_p = 6 * K_user * f_q1 * f_q2
        d_p = C0

        if force_integer_coeffs:
            a_final = a_p.limit_denominator(1)
            if a_final == 0: 
                raise ValueError("Mit diesen Extrema-Vorgaben und 'Ganzzahlige Koeffizienten' "
                                 "wird der x^3-Term Null. Kein kubisches Polynom möglich. "
                                 "Ändern Sie die y-Werte oder deaktivieren Sie die Option.")
            return cls(a_final, 
                       b_p.limit_denominator(1), 
                       c_p.limit_denominator(1), 
                       d_p.limit_denominator(1))
        else:
            if a_p == 0 : 
                 raise ValueError("Interner Fehler: a_p wurde Null, obwohl K_user nicht Null sein sollte (keine Ganzzahligkeit erzwungen).")
            return cls(a_p, b_p, c_p, d_p)


    @classmethod
    def from_inflection_point_and_delta_q(cls, xw, delta_q, K_user, C0):
        f_xw,f_delta_q,f_K_user,f_C0 = Fraction(xw),Fraction(delta_q),Fraction(K_user),Fraction(C0)
        if f_delta_q == 0: raise ValueError("delta_q muss ungleich Null sein.")
        if f_K_user == 0: raise ValueError("K_user darf nicht Null sein.")
        C1_val = 3 * f_K_user * (f_xw**2 - f_delta_q**2)
        a_p,b_p,c_p,d_p = f_K_user, -3*f_K_user*f_xw, C1_val, f_C0
        return cls(a_p, b_p, c_p, d_p)

    @classmethod
    def from_integer_inflection_real_extrema(cls, xw, D0, a_coeff, C1_coeff):
        f_xw,f_D0,f_a_coeff,f_C1_coeff = Fraction(xw),Fraction(D0),Fraction(a_coeff),Fraction(C1_coeff)
        if f_a_coeff == 0: raise ValueError("a_coeff darf nicht Null sein.")
        if f_C1_coeff != 0 and (f_C1_coeff * f_a_coeff) >= 0: 
            raise ValueError("a_coeff und C1_coeff müssen unterschiedliche Vorzeichen haben.")
        A,B,C,D = f_a_coeff, -3*f_a_coeff*f_xw, 3*f_a_coeff*f_xw**2+f_C1_coeff, -f_a_coeff*f_xw**3-f_C1_coeff*f_xw+f_D0
        return cls(A,B,C,D)

    @classmethod
    def from_zero_and_two_integer_roots(cls, r2, r3, scale_factor_S):
        f_r2,f_r3,f_S = Fraction(r2),Fraction(r3),Fraction(scale_factor_S)
        if f_S == 0: raise ValueError("scale_factor_S darf nicht Null sein.")
        a_p,b_p,c_p,d_p = f_S, -f_S*(f_r2+f_r3), f_S*f_r2*f_r3, Fraction(0)
        return cls(a_p,b_p,c_p,d_p)
        
    @classmethod
    def from_quadratic_vertex_and_scale(cls, h, k, scale_factor_sq, force_integer_coeffs=False):
        f_h, f_k, f_scale = Fraction(h), Fraction(k), Fraction(scale_factor_sq)
        if f_scale == 0: raise ValueError("scale_factor_sq darf nicht Null sein.")
        b_p, c_p, d_p = f_scale, -2*f_scale*f_h, f_scale*f_h**2+f_k
        if force_integer_coeffs:
            b_final = b_p.limit_denominator(1)
            if b_final == 0:
                raise ValueError("Mit diesen Scheitelpunkt-Vorgaben und 'Ganzzahlige Koeffizienten' "
                                 "wird der x^2-Term Null. Kein quadratisches Polynom möglich. "
                                 "Ändern Sie die Werte oder deaktivieren Sie die Option.")
            return cls(0, b_final, c_p.limit_denominator(1), d_p.limit_denominator(1))
        else:
            return cls(0, b_p, c_p, d_p)

    @classmethod
    def from_linear_root_and_y_intercept(cls, root_x, y_intercept):
        f_root_x, f_y_intercept = Fraction(root_x), Fraction(y_intercept)
        d_p = f_y_intercept
        if f_root_x == 0:
            if f_y_intercept != 0: raise ValueError("Wenn Nullstelle x=0, muss y-Abschnitt auch 0 sein.")
            else: raise ValueError("Steigung unbestimmt (0/0). Koeffizienten direkt eingeben.")
        else: c_p = -f_y_intercept / f_root_x
        return cls(0, 0, c_p, d_p)

    def _find_quadratic_roots(self, A_frac, B_frac, C_frac):
        A, B, C_coeff = float(A_frac), float(B_frac), float(C_frac)
        if abs(A) < 1e-9: 
            if abs(B) < 1e-9: return [] 
            return [-C_coeff / B]
        discriminant = B**2 - 4*A*C_coeff
        if discriminant < -1e-9: return [] 
        elif abs(discriminant) < 1e-9: return [-B / (2*A)]
        else:
            sqrt_D = math.sqrt(max(0, discriminant)) 
            return sorted([(-B + sqrt_D)/(2*A), (-B - sqrt_D)/(2*A)])

    def find_critical_points_x(self):
        if self.a == 0 : 
            if self.b == 0: return [] 
            A_prime, B_prime, C_prime = self.derivative1_coeffs() 
            return self._find_quadratic_roots(A_prime, B_prime, C_prime) 
        A_prime, B_prime, C_prime = self.derivative1_coeffs()
        return self._find_quadratic_roots(A_prime, B_prime, C_prime)

    def find_inflection_point_x(self):
        if self.a == 0: return None 
        return float(-self.b / (3 * self.a)) 

    def find_roots(self, integer_check_range=20):
        roots = []
        if self.a == 0 and self.b == 0 and self.c == 0 and self.d == 0: return ["Alle reellen Zahlen"]
        if self.a == 0 and self.b == 0 and self.c == 0 and self.d != 0: return [] 
        if self.a == 0 and self.b == 0 and self.c != 0: return [float(-self.d / self.c)]
        if self.a == 0 and self.b != 0: return self._find_quadratic_roots(self.b, self.c, self.d)
        if self.d == 0:
            roots.append(0.0)
            quad_roots = self._find_quadratic_roots(self.a, self.b, self.c)
            for r in quad_roots:
                if not any(abs(r - er) < 1e-7 for er in roots): roots.append(r)
            return sorted(list(set(roots))) 
        for i in range(1, integer_check_range + 1):
            for test_x_val in [i, -i]:
                if abs(self.evaluate(test_x_val)) < 1e-9: 
                    roots.append(float(test_x_val))
                    aq,bq,cq = self.a, self.b+self.a*Fraction(test_x_val), self.c+self.b*Fraction(test_x_val)+self.a*Fraction(test_x_val)**2
                    quad_roots = self._find_quadratic_roots(aq,bq,cq)
                    for r_q in quad_roots: 
                         if not any(abs(r_q - er) < 1e-7 for er in roots): roots.append(r_q)
                    return sorted(list(set(roots))) 
        try:
            np_coeffs = [float(c) for c in [self.a, self.b, self.c, self.d]]
            np_roots = np.roots(np_coeffs)
            real_roots = [r.real for r in np_roots if abs(r.imag) < 1e-7]
            unique_real_roots = []
            for r_np in sorted(real_roots): 
                if not any(abs(r_np - ur) < 1e-7 for ur in unique_real_roots): unique_real_roots.append(r_np)
            return unique_real_roots
        except Exception as e:
            print(f"Fehler bei numpy.roots: {e}")
            return ["Numerische Lösung fehlgeschlagen"]

class PolynomialAppGUI:
    MAX_PARAM_ROWS = 5 # Maximale Anzahl von Parameterzeilen, die wir erwarten

    def __init__(self, master):
        self.master = master
        master.title("Polynomfunktionen Generator")
        master.geometry("800x950") 

        self.style = ttk.Style()
        self.style.theme_use('clam') 

        self.paned_window = ttk.PanedWindow(master, orient=tk.VERTICAL)
        self.paned_window.pack(expand=True, fill=tk.BOTH, padx=10, pady=10)

        top_frame = ttk.Frame(self.paned_window, padding="5")
        self.paned_window.add(top_frame, weight=2) 

        plot_outer_frame = ttk.LabelFrame(self.paned_window, text="Graphische Darstellung", padding="5")
        self.paned_window.add(plot_outer_frame, weight=3) 

        poly_type_frame = ttk.LabelFrame(top_frame, text="Polynomtyp", padding="10")
        poly_type_frame.pack(pady=5, fill=tk.X)
        self.poly_type_var = tk.StringVar(value="Kubisch")
        poly_types = ["Kubisch", "Quadratisch", "Linear"]
        ttk.Label(poly_type_frame, text="Typ wählen:").pack(side=tk.LEFT, padx=5)
        self.poly_type_menu = ttk.Combobox(poly_type_frame, textvariable=self.poly_type_var,
                                           values=poly_types, state="readonly", width=20)
        self.poly_type_menu.pack(side=tk.LEFT, padx=5)
        self.poly_type_menu.bind("<<ComboboxSelected>>", self.update_ui_for_poly_type)

        self.construction_frame = ttk.LabelFrame(top_frame, text="Primäre Konstruktionseigenschaft", padding="10")
        
        self.construction_method_var = tk.StringVar()
        self.cubic_construction_options = {
            "Sattelpunkt definieren": "from_saddle_point",
            "Extrema-Lagen (x-Werte) definieren": "from_extrema_locations",
            "Extrema (x- und y-Werte) definieren": "from_extrema_locations_and_values", 
            "WP & Abstand zu Extrema definieren": "from_inflection_point_and_delta_q",
            "WP (x-Wert), Extrema reell": "from_integer_inflection_real_extrema",
            "Nullstellen (0, r2, r3) definieren": "from_zero_and_two_integer_roots"
        }
        self.quadratic_construction_options = {
            "Aus Scheitelpunkt und Skalierungsfaktor": "from_quadratic_vertex_and_scale",
            "Aus Koeffizienten (b,c,d)": "direct_coeffs_quad"
        }
        self.linear_construction_options = {
            "Aus Nullstelle und y-Achsenabschnitt": "from_linear_root_and_y_intercept",
            "Aus Koeffizienten (c,d)": "direct_coeffs_lin"
        }
        self.construction_method_label = ttk.Label(self.construction_frame, text="Konstruktion durch:")
        self.construction_method_menu = ttk.Combobox(self.construction_frame, 
                                                     textvariable=self.construction_method_var,
                                                     state="readonly", width=40)
        self.construction_method_menu.bind("<<ComboboxSelected>>", self._update_construction_params_display) 
        
        self.modifier_frame = ttk.LabelFrame(top_frame, text="Globale Modifikatoren", padding="10") 
        self.force_int_coeffs_var = tk.BooleanVar(value=True)
        self.force_int_coeffs_check = ttk.Checkbutton(self.modifier_frame, 
                                             text="Ganzzahlige Koeffizienten erzwingen (wo anwendbar)",
                                             variable=self.force_int_coeffs_var,
                                             command=self.on_modifier_toggle) 
        self.force_int_coeffs_check.pack(side=tk.LEFT, padx=5)

        self.aim_int_y_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(self.modifier_frame, 
                        text="Ganzzahlige y-Werte an speziellen Punkten anstreben", 
                        variable=self.aim_int_y_var,
                        command=self.on_modifier_toggle).pack(side=tk.LEFT, padx=5)

        self.input_params_frame = ttk.LabelFrame(top_frame, text="Parameter", padding="10")
        self.input_params_frame.pack(pady=5, fill=tk.X)
        self.input_params_frame.columnconfigure(0, weight=1) 
        self.input_params_frame.columnconfigure(1, weight=2) 
        
        # Pre-create parameter rows
        self.param_row_widgets = []
        for i in range(self.MAX_PARAM_ROWS):
            lbl = ttk.Label(self.input_params_frame, text="")
            entry_var = tk.StringVar()
            entry = ttk.Entry(self.input_params_frame, textvariable=entry_var, width=30)
            hint_lbl = ttk.Label(self.input_params_frame, text="", foreground="blue", font=("TkDefaultFont", 8))
            self.param_row_widgets.append({'label': lbl, 'entry_var': entry_var, 'entry': entry, 'hint': hint_lbl})
            # Widgets werden erst in _configure_param_entry_row gegridded und sichtbar gemacht

        self.param_entries = {} # Dictionary to map param_name to its StringVar for get_param_value

        button_frame = ttk.Frame(top_frame)
        button_frame.pack(pady=10)

        generate_button = ttk.Button(button_frame, text="Polynom generieren und analysieren", command=self.generate_and_display)
        generate_button.pack(side=tk.LEFT, padx=5)
        
        lucky_button = ttk.Button(button_frame, text="Ich fühle mich glücklich!", command=self.i_feel_lucky)
        lucky_button.pack(side=tk.LEFT, padx=5)

        output_frame = ttk.LabelFrame(top_frame, text="Ergebnisse", padding="10")
        output_frame.pack(pady=5, expand=True, fill=tk.BOTH)
        output_frame.columnconfigure(1, weight=1) 

        self.output_fields = {}
        output_labels_texts = {
            "P(x)": "Polynom P(x):", "P'(x)": "Ableitung P'(x):", "P''(x)": "Ableitung P''(x):",
            "roots": "Nullstellen P(x):",
            "crit_x": "Kritische Stellen x:", "crit_y": "Y-Werte an krit. Stellen:",
            "infl_x": "Wendestelle x:", "infl_y": "Y-Wert am Wendepunkt:",
            "saddle_info": "Sattelpunkt Info:"
        }
        row_idx = 0
        for key, text in output_labels_texts.items():
            ttk.Label(output_frame, text=text).grid(row=row_idx, column=0, sticky=tk.NW, padx=5, pady=2)
            text_widget = tk.Text(output_frame, height=1, width=60, wrap=tk.WORD, relief=tk.SUNKEN, borderwidth=1, font=("Consolas", 10))
            if key in ["P(x)", "P'(x)", "P''(x)"]: text_widget.config(height=2)
            text_widget.grid(row=row_idx, column=1, sticky=tk.EW, padx=5, pady=2)
            text_widget.config(state=tk.DISABLED) 
            self.output_fields[key] = text_widget
            row_idx += 1
        
        self.fig = Figure(figsize=(7, 5), dpi=100) 
        self.plot_ax = self.fig.add_subplot(111)
        self.plot_ax.set_xlabel("x")
        self.plot_ax.set_ylabel("P(x)")
        self.plot_ax.grid(True, linestyle=':', alpha=0.7)
        self.plot_ax.axhline(0, color='black', linewidth=0.5)
        self.plot_ax.axvline(0, color='black', linewidth=0.5)

        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_outer_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.canvas.draw() 
            
        self.update_ui_for_poly_type() 

    def on_modifier_toggle(self):
        self._apply_param_hints()

    def update_ui_for_poly_type(self, event=None):
        poly_type = self.poly_type_var.get()
        self.clear_param_input_fields_display() # Nur Widgets ausblenden und Texte löschen

        if poly_type == "Kubisch":
            self.construction_frame.config(text="Primäre Konstruktionseigenschaft (Kubisch)")
            self.construction_method_label.pack(side=tk.LEFT, padx=5)
            self.construction_method_menu.config(values=list(self.cubic_construction_options.keys()))
            self.construction_method_menu.pack(side=tk.LEFT, padx=5, expand=True, fill=tk.X)
            self.construction_frame.pack(pady=5, fill=tk.X, before=self.input_params_frame) 
            self.modifier_frame.pack(pady=5, fill=tk.X, before=self.input_params_frame) 
            if not self.construction_method_var.get() or self.construction_method_var.get() not in self.cubic_construction_options:
                 self.construction_method_menu.current(0) 
            self._update_cubic_construction_params_create_fields() 
        elif poly_type == "Quadratisch":
            self.construction_frame.config(text="Konstruktionsmethode (Quadratisch)")
            self.construction_method_label.pack(side=tk.LEFT, padx=5)
            self.construction_method_menu.config(values=list(self.quadratic_construction_options.keys()))
            self.construction_method_menu.pack(side=tk.LEFT, padx=5, expand=True, fill=tk.X)
            self.construction_frame.pack(pady=5, fill=tk.X, before=self.input_params_frame)
            self.modifier_frame.pack(pady=5, fill=tk.X, before=self.input_params_frame) 
            if not self.construction_method_var.get() or self.construction_method_var.get() not in self.quadratic_construction_options:
                 self.construction_method_menu.current(0)
            self._update_quadratic_construction_params_create_fields() 
        else: # Linear
            self.construction_frame.config(text="Konstruktionsmethode (Linear)")
            self.construction_method_label.pack(side=tk.LEFT, padx=5)
            self.construction_method_menu.config(values=list(self.linear_construction_options.keys()))
            self.construction_method_menu.pack(side=tk.LEFT, padx=5, expand=True, fill=tk.X)
            self.construction_frame.pack(pady=5, fill=tk.X, before=self.input_params_frame)
            self.modifier_frame.pack_forget() 
            if not self.construction_method_var.get() or self.construction_method_var.get() not in self.linear_construction_options:
                 self.construction_method_menu.current(0)
            self._update_linear_construction_params_create_fields() 


    def clear_param_input_fields_display(self):
        """Blendet alle vorab erstellten Parameterzeilen aus und löscht ihre Werte/Texte."""
        for i in range(self.MAX_PARAM_ROWS):
            row = self.param_row_widgets[i]
            row['label'].grid_remove()
            row['entry'].grid_remove()
            row['hint'].grid_remove()
            row['label'].config(text="")
            row['entry_var'].set("")
            row['hint'].config(text="")
        self.param_entries.clear() # Das Dictionary, das param_name auf entry_var mappt

    def _configure_param_entry_row(self, row_index, param_name, label_text, default_value=""):
        """Konfiguriert und zeigt eine vorab erstellte Parameterzeile an."""
        if row_index >= self.MAX_PARAM_ROWS:
            print(f"Warnung: Es wird versucht, mehr als {self.MAX_PARAM_ROWS} Parameterzeilen zu verwenden.")
            return row_index # Nichts tun, wenn Index außerhalb des Bereichs liegt

        row = self.param_row_widgets[row_index]
        
        row['label'].config(text=label_text)
        row['entry_var'].set(str(default_value))
        # Hinweis wird später durch _apply_param_hints gesetzt

        row['label'].grid(row=row_index, column=0, sticky=tk.W, padx=5, pady=3)
        row['entry'].grid(row=row_index, column=1, sticky=tk.EW, padx=5, pady=3)
        row['hint'].grid(row=row_index, column=2, sticky=tk.W, padx=5)
        
        self.param_entries[param_name] = row['entry_var'] # Wichtig für get_param_value
        return row_index + 1
    
    def _update_construction_params_display(self, event=None):
        poly_type = self.poly_type_var.get()
        if poly_type == "Kubisch":
            self._update_cubic_construction_params_create_fields()
        elif poly_type == "Quadratisch":
            self._update_quadratic_construction_params_create_fields()
        elif poly_type == "Linear":
            self._update_linear_construction_params_create_fields()

    def _apply_param_hints(self):
        poly_type = self.poly_type_var.get()
        construction_key = self.construction_method_var.get()
        force_int_active = self.force_int_coeffs_var.get()
        aim_int_y_active = self.aim_int_y_var.get()

        # Zuerst alle (sichtbaren) Hints leeren
        for i in range(self.MAX_PARAM_ROWS):
            row = self.param_row_widgets[i]
            if row['hint'].winfo_ismapped(): # Nur wenn der Hint sichtbar ist
                 row['hint'].config(text="")


        if poly_type == "Kubisch":
            method_name = self.cubic_construction_options.get(construction_key)
            if method_name == "from_saddle_point":
                self.param_row_widgets[2]['hint'].config(text="Ganze Zahl & Vielfaches von 6" if force_int_active else "Kann Bruch sein (z.B. 1/6)")
                if aim_int_y_active: self.param_row_widgets[1]['hint'].config(text="Ganzzahlig für ganzz. y am SP")
            elif method_name == "from_extrema_locations":
                if aim_int_y_active: self.param_row_widgets[3]['hint'].config(text="Ganzzahlig für ganzz. y-Extrema") # C0 ist 4. Parameter (Index 3)
            elif method_name == "from_extrema_locations_and_values":
                if aim_int_y_active: 
                    self.param_row_widgets[1]['hint'].config(text="Ganzzahliger y-Wert") # y1_val
                    self.param_row_widgets[3]['hint'].config(text="Ganzzahliger y-Wert") # y2_val
                if force_int_active: 
                     self.param_row_widgets[0]['hint'].config(text="Ganzz. Koeff. werden angestrebt") # q1_val
            elif method_name == "from_inflection_point_and_delta_q":
                if aim_int_y_active: self.param_row_widgets[3]['hint'].config(text="Ganzzahlig für ganzz. y-Werte (beeinflusst d)")
            elif method_name == "from_integer_inflection_real_extrema":
                self.param_row_widgets[3]['hint'].config(text="Muss anderes Vorzeichen als 'a' haben für Extrema") # C1_coeff
                if aim_int_y_active: self.param_row_widgets[1]['hint'].config(text="Ganzzahlig für ganzz. y am WP") # D0
        
        elif poly_type == "Quadratisch":
            method_name = self.quadratic_construction_options.get(construction_key)
            if method_name == "from_quadratic_vertex_and_scale":
                if aim_int_y_active: self.param_row_widgets[1]['hint'].config(text="Ganzzahlig für ganzz. y am Scheitel") # k_vertex
                if force_int_active: self.param_row_widgets[2]['hint'].config(text="Ganzz. Koeff. wenn h,k,scale_sq ganzz. oder passende Brüche") #scale_sq
        
        elif poly_type == "Linear":
            method_name = self.linear_construction_options.get(construction_key)
            if method_name == "from_linear_root_and_y_intercept":
                 self.param_row_widgets[0]['hint'].config(text="Wenn 0, muss y-Abschnitt 0 sein (Steigung dann unbestimmt)") # root_x_lin

    def _update_cubic_construction_params_create_fields(self):
        self.clear_param_input_fields_display() 
        self.input_params_frame.config(text="Parameter für Kubische Konstruktion") 
        construction_key = self.construction_method_var.get()
        method_name = self.cubic_construction_options.get(construction_key)
        current_row = 0
        if method_name == "from_saddle_point":
            current_row = self._configure_param_entry_row(current_row, "x0", "x-Sattelpunkt (x0):", "1")
            current_row = self._configure_param_entry_row(current_row, "D0", "y-Sattelpunkt (D0):", "2")
            current_row = self._configure_param_entry_row(current_row, "P_triple_prime_0", "P'''(x0):", "6")
        elif method_name == "from_extrema_locations":
            current_row = self._configure_param_entry_row(current_row, "q1", "x-Extremum 1 (q1):", "0")
            current_row = self._configure_param_entry_row(current_row, "q2", "x-Extremum 2 (q2):", "2")
            current_row = self._configure_param_entry_row(current_row, "K_user", "Skalierungsfaktor K:", "1")
            current_row = self._configure_param_entry_row(current_row, "C0", "Konstante C0 (y-Versch.):", "0")
        elif method_name == "from_extrema_locations_and_values": 
            current_row = self._configure_param_entry_row(current_row, "q1_val", "x-Extremum 1 (q1):", "0")
            current_row = self._configure_param_entry_row(current_row, "y1_val", "y-Wert bei q1 (y1):", "0")
            current_row = self._configure_param_entry_row(current_row, "q2_val", "x-Extremum 2 (q2):", "2")
            current_row = self._configure_param_entry_row(current_row, "y2_val", "y-Wert bei q2 (y2):", "4")
        elif method_name == "from_inflection_point_and_delta_q":
            current_row = self._configure_param_entry_row(current_row, "xw", "x-Wendepunkt (xw):", "1")
            current_row = self._configure_param_entry_row(current_row, "delta_q", "Abstand zu Extrema (Δq):", "1")
            current_row = self._configure_param_entry_row(current_row, "K_user", "Skalierungsfaktor K:", "1")
            current_row = self._configure_param_entry_row(current_row, "C0", "Konstante C0 (y-Versch.):", "0")
        elif method_name == "from_integer_inflection_real_extrema":
            current_row = self._configure_param_entry_row(current_row, "xw", "x-Wendepunkt (xw):", "0")
            current_row = self._configure_param_entry_row(current_row, "D0", "y-Wendepunkt (D0):", "0")
            current_row = self._configure_param_entry_row(current_row, "a_coeff", "Koeff. a (von x^3):", "1")
            current_row = self._configure_param_entry_row(current_row, "C1_coeff", "Koeff. C1 (linearer Term um WP):", "-3")
        elif method_name == "from_zero_and_two_integer_roots":
            current_row = self._configure_param_entry_row(current_row, "r2", "Nullstelle r2:", "1")
            current_row = self._configure_param_entry_row(current_row, "r3", "Nullstelle r3:", "2")
            current_row = self._configure_param_entry_row(current_row, "scale_factor_S", "Skalierungsfaktor S:", "1")
        
        self._apply_param_hints()


    def _update_quadratic_construction_params_create_fields(self):
        self.clear_param_input_fields_display()
        self.input_params_frame.config(text="Parameter für Quadratische Konstruktion")
        construction_key = self.construction_method_var.get()
        method_name = self.quadratic_construction_options.get(construction_key)
        current_row = 0
        if method_name == "from_quadratic_vertex_and_scale":
            current_row = self._configure_param_entry_row(current_row, "h_vertex", "Scheitelpunkt x (h):", "1")
            current_row = self._configure_param_entry_row(current_row, "k_vertex", "Scheitelpunkt y (k):", "2")
            current_row = self._configure_param_entry_row(current_row, "scale_sq", "Skalierungsfaktor:", "1")
        elif method_name == "direct_coeffs_quad":
            current_row = self._configure_param_entry_row(current_row, "b_quad", "Koeffizient b (für x^2):", "1")
            current_row = self._configure_param_entry_row(current_row, "c_quad", "Koeffizient c (für x):", "0")
            current_row = self._configure_param_entry_row(current_row, "d_quad", "Koeffizient d (Konstante):", "0")
        
        self._apply_param_hints()


    def _update_linear_construction_params_create_fields(self):
        self.clear_param_input_fields_display()
        self.input_params_frame.config(text="Parameter für Lineare Konstruktion")
        construction_key = self.construction_method_var.get()
        method_name = self.linear_construction_options.get(construction_key)
        current_row = 0
        if method_name == "from_linear_root_and_y_intercept":
            current_row = self._configure_param_entry_row(current_row, "root_x_lin", "Nullstelle x:", "1")
            current_row = self._configure_param_entry_row(current_row, "y_intercept_lin", "y-Achsenabschnitt (d):", "2")
        elif method_name == "direct_coeffs_lin":
            current_row = self._configure_param_entry_row(current_row, "c_lin", "Koeffizient c (Steigung m):", "1")
            current_row = self._configure_param_entry_row(current_row, "d_lin", "Koeffizient d (y-Abschnitt n):", "0")
        
        self._apply_param_hints()


    def get_param_value(self, param_name, param_type_hint=Fraction):
        try:
            # param_entries enthält jetzt die StringVars der sichtbaren/konfigurierten Felder
            if param_name not in self.param_entries:
                raise ValueError(f"Parameter '{param_name}' ist für die aktuelle Auswahl nicht verfügbar.")
            val_str = self.param_entries[param_name].get().strip()
            if not val_str: raise ValueError(f"Parameter '{param_name}' darf nicht leer sein.")
            
            f_val = Fraction(val_str)
            if param_type_hint == int:
                if f_val.denominator != 1: raise ValueError(f"'{param_name}' muss ganze Zahl sein, nicht '{val_str}'.")
                return f_val.numerator
            return f_val 
        except KeyError: 
            # Dieser Fall sollte durch die obige Prüfung abgedeckt sein, aber als Fallback
            raise ValueError(f"Parameter '{param_name}' nicht im UI gefunden (Programmierfehler).")
        except (ValueError, TypeError) as e: 
            raise ValueError(f"Ungültiger Wert für '{param_name}': '{val_str}' ({e}). Erwartet Zahl oder Bruch (z.B. 1/3).")


    def generate_and_display(self):
        for key in self.output_fields:
            self.output_fields[key].config(state=tk.NORMAL); self.output_fields[key].delete(1.0, tk.END); self.output_fields[key].config(state=tk.DISABLED)
        try:
            poly_type = self.poly_type_var.get()
            poly = None
            force_int_coeffs_global = self.force_int_coeffs_var.get()

            if poly_type == "Kubisch":
                construction_key, method_name = self.construction_method_var.get(), self.cubic_construction_options.get(self.construction_method_var.get())
                if not construction_key: messagebox.showerror("Auswahlfehler", "Bitte Konstruktionsmethode wählen."); return
                
                if method_name == "from_saddle_point":
                    x0,D0,P_ppp = self.get_param_value("x0"), self.get_param_value("D0"), self.get_param_value("P_triple_prime_0", int if force_int_coeffs_global else Fraction)
                    poly = CubicPolynomial.from_saddle_point(x0,D0,P_ppp,force_int_coeffs_global)
                elif method_name == "from_extrema_locations":
                    q1,q2,K,C0 = self.get_param_value("q1"),self.get_param_value("q2"),self.get_param_value("K_user"),self.get_param_value("C0")
                    poly = CubicPolynomial.from_extrema_locations(q1,q2,K,C0)
                elif method_name == "from_extrema_locations_and_values": 
                    q1,y1,q2,y2 = self.get_param_value("q1_val"),self.get_param_value("y1_val"),self.get_param_value("q2_val"),self.get_param_value("y2_val")
                    poly = CubicPolynomial.from_extrema_locations_and_values(q1,q2,y1,y2,force_int_coeffs_global)
                elif method_name == "from_inflection_point_and_delta_q":
                    xw,dq,K,C0 = self.get_param_value("xw"),self.get_param_value("delta_q"),self.get_param_value("K_user"),self.get_param_value("C0")
                    poly = CubicPolynomial.from_inflection_point_and_delta_q(xw,dq,K,C0)
                elif method_name == "from_integer_inflection_real_extrema":
                    xw,D0,a_c,C1_c = self.get_param_value("xw"),self.get_param_value("D0"),self.get_param_value("a_coeff"),self.get_param_value("C1_coeff")
                    poly = CubicPolynomial.from_integer_inflection_real_extrema(xw,D0,a_c,C1_c)
                elif method_name == "from_zero_and_two_integer_roots":
                    r2,r3,S = self.get_param_value("r2"),self.get_param_value("r3"),self.get_param_value("scale_factor_S")
                    poly = CubicPolynomial.from_zero_and_two_integer_roots(r2,r3,S)
            
            elif poly_type == "Quadratisch":
                construction_key_quad, method_name_quad = self.construction_method_var.get(), self.quadratic_construction_options.get(self.construction_method_var.get())
                if not construction_key_quad: messagebox.showerror("Auswahlfehler", "Bitte Konstruktionsmethode wählen."); return
                if method_name_quad == "from_quadratic_vertex_and_scale":
                    h,k,scale_sq = self.get_param_value("h_vertex"),self.get_param_value("k_vertex"),self.get_param_value("scale_sq")
                    poly = CubicPolynomial.from_quadratic_vertex_and_scale(h,k,scale_sq,force_int_coeffs_global)
                elif method_name_quad == "direct_coeffs_quad":
                    b,c,d = self.get_param_value("b_quad"),self.get_param_value("c_quad"),self.get_param_value("d_quad")
                    poly = CubicPolynomial(0,b,c,d) 
            
            elif poly_type == "Linear":
                construction_key_lin, method_name_lin = self.construction_method_var.get(), self.linear_construction_options.get(self.construction_method_var.get())
                if not construction_key_lin: messagebox.showerror("Auswahlfehler", "Bitte Konstruktionsmethode wählen."); return
                if method_name_lin == "from_linear_root_and_y_intercept":
                    root_x,y_int = self.get_param_value("root_x_lin"),self.get_param_value("y_intercept_lin")
                    poly = CubicPolynomial.from_linear_root_and_y_intercept(root_x,y_int)
                elif method_name_lin == "direct_coeffs_lin":
                    c,d = self.get_param_value("c_lin"),self.get_param_value("d_lin")
                    poly = CubicPolynomial(0,0,c,d) 

            if poly: self.display_poly_info(poly)
            else: messagebox.showerror("Fehler", "Polynom konnte nicht generiert werden.")
        except ValueError as e: messagebox.showerror("Eingabefehler", str(e))
        except Exception as e: messagebox.showerror("Unerwarteter Fehler", f"Fehler: {str(e)}")

    def display_poly_info(self, poly):
        def format_deriv_poly_str(coeffs_frac):
            if len(coeffs_frac) == 3: temp_poly = CubicPolynomial(0, coeffs_frac[0], coeffs_frac[1], coeffs_frac[2])
            elif len(coeffs_frac) == 2: temp_poly = CubicPolynomial(0, 0, coeffs_frac[0], coeffs_frac[1])
            elif len(coeffs_frac) == 1: temp_poly = CubicPolynomial(0,0,0,coeffs_frac[0])
            else: return "N/A"
            return str(temp_poly)
        p_prime_coeffs_frac, p_double_prime_coeffs_frac = poly.derivative1_coeffs(), poly.derivative2_coeffs()
        roots_found = poly.find_roots()
        roots_str = (roots_found[0] if isinstance(roots_found,list) and roots_found and isinstance(roots_found[0],str) and "Alle" in roots_found[0] 
                     else (roots_found[0] if isinstance(roots_found,list) and roots_found and isinstance(roots_found[0],str)
                           else (", ".join([f"{float(r):.3f}" for r in roots_found]) if roots_found else "Keine (reellen) gefunden")))
        crit_x_float = poly.find_critical_points_x()
        crit_y_str = ", ".join([f"{float(poly.evaluate(x)):.3f}" for x in crit_x_float]) if crit_x_float else "Keine"
        crit_x_str = ", ".join([f"{x:.3f}" for x in crit_x_float]) if crit_x_float else "Keine"
        infl_x_float = poly.find_inflection_point_x()
        infl_y_str, infl_x_str, saddle_info_str = "N/A", "N/A", "Kein Sattelpunkt"
        if infl_x_float is not None:
            infl_x_str, infl_y_str = f"{infl_x_float:.3f}", f"{float(poly.evaluate(infl_x_float)):.3f}"
            if crit_x_float and len(crit_x_float) == 1 and abs(crit_x_float[0] - infl_x_float) < 1e-7:
                saddle_info_str = f"Ja, bei x = {infl_x_float:.3f}"
        outputs = {"P(x)": str(poly), "P'(x)": format_deriv_poly_str(p_prime_coeffs_frac),
                   "P''(x)": format_deriv_poly_str(p_double_prime_coeffs_frac), "roots": roots_str,
                   "crit_x": crit_x_str, "crit_y": crit_y_str, "infl_x": infl_x_str,
                   "infl_y": infl_y_str, "saddle_info": saddle_info_str}
        for key, value in outputs.items():
            self.output_fields[key].config(state=tk.NORMAL); self.output_fields[key].delete(1.0, tk.END)
            self.output_fields[key].insert(tk.END, value); self.output_fields[key].config(state=tk.DISABLED)

        self.plot_ax.clear() 
        x_min_plot, x_max_plot = -5, 5 
        all_important_x_plot = [float(r) for r in roots_found if not isinstance(r, str)] + \
                               [float(x) for x in crit_x_float]
        if infl_x_float is not None: all_important_x_plot.append(infl_x_float)
        
        if all_important_x_plot:
            min_x_val, max_x_val = min(all_important_x_plot), max(all_important_x_plot) 
            padding = (max_x_val - min_x_val) * 0.2 if (max_x_val - min_x_val) > 1e-6 else 2.0
            x_min_plot, x_max_plot = min_x_val - padding, max_x_val + padding
        
        x_vals = np.linspace(x_min_plot, x_max_plot, 400)
        y_vals = poly.evaluate(x_vals) 

        self.plot_ax.plot(x_vals, np.array(y_vals, dtype=float), label=f"P(x) = {str(poly)}", color='blue')
        
        if roots_found and not isinstance(roots_found[0], str):
            self.plot_ax.scatter([float(r) for r in roots_found], [0]*len(roots_found), color='black', s=50, zorder=5, label="Nullstellen")

        if crit_x_float:
            is_saddle_for_plot = infl_x_float is not None and len(crit_x_float) == 1 and abs(crit_x_float[0] - infl_x_float) < 1e-7
            non_saddle_crit_x = [x for x in crit_x_float if not (is_saddle_for_plot and abs(x - infl_x_float) < 1e-7)]
            non_saddle_crit_y = [float(poly.evaluate(x)) for x in non_saddle_crit_x]
            if non_saddle_crit_x:
                 self.plot_ax.scatter(non_saddle_crit_x, non_saddle_crit_y, color='red', zorder=5, s=70, label="Extrema", edgecolors='black')
        
        if infl_x_float is not None:
            is_saddle_for_plot = len(crit_x_float) == 1 and abs(crit_x_float[0] - infl_x_float) < 1e-7
            color_infl = 'purple' if is_saddle_for_plot else 'green'
            label_infl = "Sattelpunkt" if is_saddle_for_plot else "Wendepunkt"
            marker_infl = '*' if is_saddle_for_plot else 'P' 
            size_infl = 200 if is_saddle_for_plot else 150
            self.plot_ax.scatter([infl_x_float], [float(poly.evaluate(infl_x_float))], color=color_infl, zorder=6, s=size_infl, marker=marker_infl, label=label_infl, edgecolors='black', linewidth=1)

        self.plot_ax.set_title(f"Graph von P(x)")
        self.plot_ax.set_xlabel("x")
        self.plot_ax.set_ylabel("P(x)")
        self.plot_ax.grid(True, linestyle=':', alpha=0.7)
        self.plot_ax.axhline(0, color='black', linewidth=0.5)
        self.plot_ax.axvline(0, color='black', linewidth=0.5)
        
        if y_vals.size > 0:
            y_vals_float = np.array(y_vals, dtype=float)
            min_y_plot, max_y_plot = np.min(y_vals_float), np.max(y_vals_float)
            y_range_plot = max_y_plot - min_y_plot
            padding_y = y_range_plot * 0.1 if abs(y_range_plot) > 1e-6 else 1.0
            self.plot_ax.set_ylim(min_y_plot - padding_y, max_y_plot + padding_y)

        self.plot_ax.legend(loc='best')
        self.canvas.draw()
        
    def _generate_random_value(self, is_integer_only=False, non_zero=False, 
                               min_val=-5, max_val=5, 
                               fraction_den_min=1, fraction_den_max=3, force_different_from=None):
        """Generiert einen zufälligen Wert (ganze Zahl oder Bruch als String)."""
        attempt_count = 0
        while attempt_count < 20: 
            if is_integer_only or random.choice([True, False]): 
                value = random.randint(min_val, max_val)
                if non_zero and value == 0:
                    attempt_count += 1; continue
                if force_different_from is not None and value == force_different_from:
                    attempt_count += 1; continue
                return str(value)
            else: 
                num_min_actual = min_val
                num_max_actual = max_val

                if non_zero:
                    if min_val == 0 and max_val == 0: 
                        num = random.choice([-1, 1]) 
                    elif min_val == 0:
                        num_min_actual = 1
                    elif max_val == 0:
                        num_max_actual = -1
                
                num = random.randint(num_min_actual, num_max_actual)
                if non_zero and num == 0: 
                    possible_nums = [i for i in range(num_min_actual, num_max_actual + 1) if i != 0]
                    if not possible_nums: 
                        num = 1 if num_max_actual >=1 else -1 
                    else:
                        num = random.choice(possible_nums)
                
                den = random.randint(fraction_den_min, fraction_den_max)
                if den == 0 : den = 1 

                val_str = f"{num}/{den}"
                try:
                    current_frac = Fraction(val_str)
                    if force_different_from is not None and current_frac == Fraction(force_different_from):
                        attempt_count +=1; continue
                    return val_str
                except ZeroDivisionError: 
                    attempt_count +=1; continue
        fallback_val = random.randint(min_val,max_val)
        while non_zero and fallback_val == 0: fallback_val = random.randint(min_val,max_val)
        while force_different_from is not None and fallback_val == force_different_from: fallback_val = random.randint(min_val,max_val)
        return str(fallback_val) 


    def i_feel_lucky(self):
        poly_type = self.poly_type_var.get()
        
        current_construction_options = {}
        if poly_type == "Kubisch":
            current_construction_options = self.cubic_construction_options
        elif poly_type == "Quadratisch":
            current_construction_options = self.quadratic_construction_options
        elif poly_type == "Linear":
            current_construction_options = self.linear_construction_options
        
        if not current_construction_options:
            self.generate_and_display() 
            return

        max_retries = 10
        for attempt in range(max_retries):
            try:
                chosen_method_key = random.choice(list(current_construction_options.keys()))
                self.construction_method_var.set(chosen_method_key) 
                
                self._update_construction_params_display() 

                self.force_int_coeffs_var.set(random.choice([True, False]))
                self.aim_int_y_var.set(random.choice([True, False]))
                self._apply_param_hints() # Nur Hints aktualisieren, nicht Felder neu erstellen
                
                x_coord_range = (-3, 3)
                y_coord_range = (-5, 5)
                scale_val_range = (-3, 3) 
                delta_range = (1, 2)
                p_ppp_forced_int_choices = [-18, -12, -6, 6, 12, 18]

                for param_name, entry_var_tk in self.param_entries.items():
                    random_val_str = "1" 
                    is_int_only = False
                    non_zero_needed = False
                    current_min, current_max = x_coord_range 

                    if param_name in ["x0", "q1", "q2", "q1_val", "q2_val", "xw", "h_vertex", "root_x_lin", "r2", "r3"]:
                        is_int_only = True 
                    elif param_name in ["D0", "C0", "k_vertex", "y_intercept_lin", "y1_val", "y2_val"]:
                         is_int_only = self.aim_int_y_var.get() 
                         current_min, current_max = y_coord_range
                    elif param_name in ["K_user", "a_coeff", "scale_factor_S", "scale_sq", "c_lin"]:
                        non_zero_needed = True
                        current_min, current_max = scale_val_range
                    elif param_name == "P_triple_prime_0":
                        if self.force_int_coeffs_var.get():
                            random_val_str = str(random.choice(p_ppp_forced_int_choices))
                            entry_var_tk.set(random_val_str)
                            continue 
                        else: 
                            non_zero_needed = True
                            current_min, current_max = scale_val_range
                    elif param_name == "delta_q":
                        is_int_only = True
                        non_zero_needed = True
                        current_min, current_max = delta_range
                    elif param_name == "C1_coeff": 
                         current_min, current_max = scale_val_range
                         
                    if param_name != "P_triple_prime_0" or not self.force_int_coeffs_var.get(): 
                        random_val_str = self._generate_random_value(is_int_only, non_zero_needed, current_min, current_max)
                    
                    entry_var_tk.set(random_val_str)

                method_key_lucky = self.construction_method_var.get() 
                
                if poly_type == "Kubisch":
                    active_method_name = self.cubic_construction_options.get(method_key_lucky)
                    if active_method_name == "from_extrema_locations" or active_method_name == "from_extrema_locations_and_values":
                        q1_var_name = "q1" if active_method_name == "from_extrema_locations" else "q1_val"
                        q2_var_name = "q2" if active_method_name == "from_extrema_locations" else "q2_val"
                        q1_var = self.param_entries.get(q1_var_name)
                        q2_var = self.param_entries.get(q2_var_name)
                        if q1_var and q2_var and q1_var.get() == q2_var.get():
                            try:
                                q1_val_int = int(Fraction(q1_var.get()))
                                q2_var.set(str(self._generate_random_value(True, False, x_coord_range[0], x_coord_range[1], force_different_from=q1_val_int)))
                            except ValueError: 
                                q2_var.set(str(self._generate_random_value(True, False, x_coord_range[0], x_coord_range[1])))
                    elif active_method_name == "from_integer_inflection_real_extrema":
                        a_coeff_var = self.param_entries.get("a_coeff")
                        c1_coeff_var = self.param_entries.get("C1_coeff")
                        if a_coeff_var and c1_coeff_var:
                            try:
                                a_val = Fraction(a_coeff_var.get())
                                c1_val = Fraction(c1_coeff_var.get())
                                if c1_val != 0 and (a_val * c1_val >= 0): 
                                    new_c1_num = -1 * random.randint(1, abs(scale_val_range[1])) if a_val > 0 else random.randint(1, abs(scale_val_range[1]))
                                    if new_c1_num == 0: new_c1_num = -1 if a_val > 0 else 1 
                                    c1_coeff_var.set(str(new_c1_num))
                                elif c1_val == 0 and a_val == 0: 
                                     a_coeff_var.set(self._generate_random_value(False, True, scale_val_range[0],scale_val_range[1]))
                            except ValueError: pass 
                elif poly_type == "Linear":
                     active_method_name = self.linear_construction_options.get(method_key_lucky)
                     if active_method_name == "from_linear_root_and_y_intercept":
                         root_x_var = self.param_entries.get("root_x_lin")
                         y_int_var = self.param_entries.get("y_intercept_lin")
                         if root_x_var and y_int_var:
                             try:
                                 if Fraction(root_x_var.get()) == 0 and Fraction(y_int_var.get()) != 0:
                                     y_int_var.set("0") 
                                 elif Fraction(root_x_var.get()) == 0 and Fraction(y_int_var.get()) == 0:
                                     root_x_var.set(self._generate_random_value(True, True, x_coord_range[0], x_coord_range[1]))
                             except ValueError: pass

                self._apply_param_hints() 
                self.generate_and_display() 
                return 
            except ValueError as e:
                print(f"I Feel Lucky - Versuch {attempt+1}/{max_retries} fehlgeschlagen: {e}")
                if attempt == max_retries - 1:
                    messagebox.showwarning("Zufallsgenerator", 
                                           f"Konnte nach {max_retries} Versuchen kein gültiges Zufallspolynom erzeugen. "
                                           "Bitte Parameter manuell anpassen oder erneut versuchen.")
        
    def _generate_random_value(self, is_integer_only=False, non_zero=False, 
                               min_val=-5, max_val=5, 
                               fraction_den_min=1, fraction_den_max=3, force_different_from=None):
        attempt_count = 0
        while attempt_count < 20: 
            if is_integer_only or random.choice([True, False]): 
                value = random.randint(min_val, max_val)
                if non_zero and value == 0:
                    attempt_count += 1; continue
                if force_different_from is not None and value == force_different_from:
                    attempt_count += 1; continue
                return str(value)
            else: 
                num_min_actual = min_val
                num_max_actual = max_val

                if non_zero:
                    if min_val == 0 and max_val == 0: 
                        num = random.choice([-1, 1]) 
                    elif min_val == 0:
                        num_min_actual = 1
                    elif max_val == 0:
                        num_max_actual = -1
                
                num = random.randint(num_min_actual, num_max_actual)
                if non_zero and num == 0: 
                    possible_nums = [i for i in range(num_min_actual, num_max_actual + 1) if i != 0]
                    if not possible_nums: 
                        num = 1 if num_max_actual >=1 else -1 
                    else:
                        num = random.choice(possible_nums)
                
                den = random.randint(fraction_den_min, fraction_den_max)
                if den == 0 : den = 1 

                val_str = f"{num}/{den}"
                try:
                    current_frac = Fraction(val_str)
                    if force_different_from is not None and current_frac == Fraction(force_different_from):
                        attempt_count +=1; continue
                    return val_str
                except ZeroDivisionError: 
                    attempt_count +=1; continue
        fallback_val = random.randint(min_val,max_val)
        while non_zero and fallback_val == 0: fallback_val = random.randint(min_val,max_val)
        while force_different_from is not None and fallback_val == force_different_from: fallback_val = random.randint(min_val,max_val)
        return str(fallback_val) 


if __name__ == '__main__':
    root = tk.Tk()
    app = PolynomialAppGUI(root)
    root.mainloop()
