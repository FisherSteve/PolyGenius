import math
import tkinter as tk
from tkinter import ttk, messagebox
from fractions import Fraction
import numpy as np # Für die allgemeine Nullstellensuche und linspace
import random # Für den "Ich fühle mich glücklich"-Button
import re # Für das Parsen von Polynom-Strings
import argparse # Für Kommandozeilenargumente

# Matplotlib-Integration für Tkinter
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

class Polynomial:
    """
    Repräsentiert ein Polynom P(x) = ex^4 + ax^3 + bx^2 + cx + d.
    Koeffizienten werden intern als Fraction-Objekte gespeichert.
    """
    def __init__(self, e, a, b, c, d): # e ist für x^4
        try:
            self.e = Fraction(e) # Koeffizient für x^4
            self.a = Fraction(a) # Koeffizient für x^3
            self.b = Fraction(b) # Koeffizient für x^2
            self.c = Fraction(c) # Koeffizient für x
            self.d = Fraction(d) # Konstanter Term
        except (TypeError, ValueError) as err_val:
            raise ValueError(f"Ungültiger Koeffiziententyp oder -wert für Polynomial: {err_val}")

    def __str__(self):
        if self.e == 0 and self.a == 0 and self.b == 0 and self.c == 0 and self.d == 0:
            return "0"
        
        processed_terms = []
        is_first_displayed_term = True
        
        coeff_power_pairs = [
            (self.e, "x^4"),
            (self.a, "x^3"),
            (self.b, "x^2"),
            (self.c, "x"),
            (self.d, "") 
        ]

        for coeff_frac, power_str in coeff_power_pairs:
            if coeff_frac == 0:
                # Überspringe den Term, außer es ist der konstante Term und alle anderen sind 0
                if not (power_str == "" and self.e == 0 and self.a == 0 and self.b == 0 and self.c == 0):
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

            # Spezialfall für konstanten Term 0, wenn es der einzige Term ist
            if power_str == "" and coeff_frac == 0 and not processed_terms and is_first_displayed_term :
                 processed_terms.append("0")
                 is_first_displayed_term = False # Verhindere, dass andere Nullen hinzugefügt werden
                 continue


            if is_first_displayed_term and val_str == "" and power_str != "": 
                processed_terms.append(f"{sign_str}{power_str}")
            else:
                processed_terms.append(f"{sign_str}{val_str}{power_str}")
            
            if coeff_frac != 0 or (power_str == "" and self.e == 0 and self.a == 0 and self.b == 0 and self.c == 0): 
                is_first_displayed_term = False
                
        final_str = "".join(processed_terms)
        if not final_str: return "0" # Fallback
        # Entferne führendes " + " falls vorhanden
        if final_str.startswith(" + "):
            final_str = final_str[3:]
        return final_str

    def evaluate(self, x): 
        # Evaluiert P(x) direkt, da die Verschiebung bereits in den Koeffizienten ist.
        if isinstance(x, np.ndarray): 
            e_f, a_f, b_f, c_f, d_f = float(self.e), float(self.a), float(self.b), float(self.c), float(self.d)
            x_float = np.array(x, dtype=float)
            return e_f * x_float**4 + a_f * x_float**3 + b_f * x_float**2 + c_f * x_float + d_f
        else: 
            x_frac = Fraction(x)
            return self.e * x_frac**4 + self.a * x_frac**3 + self.b * x_frac**2 + self.c * x_frac + self.d

    def derivative1_coeffs(self):
        # P'(x) = 4ex^3 + 3ax^2 + 2bx + c
        # Gibt Koeffizienten für die Polynomial-Klasse zurück: [e', a', b', c', d']
        return [0, 4 * self.e, 3 * self.a, 2 * self.b, self.c] 

    def derivative2_coeffs(self):
        # P''(x) = 12ex^2 + 6ax + 2b
        return [0, 0, 12 * self.e, 6 * self.a, 2 * self.b] 

    def derivative3_coeffs(self):
        # P'''(x) = 24ex + 6a
        return [0,0,0, 24 * self.e, 6 * self.a]

    @classmethod
    def from_coeffs(cls, coeffs_list):
        """ Erstellt ein Polynomial Objekt aus einer Koeffizientenliste [e, a, b, c, d] """
        e, a, b, c, d_ = Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(0)
        # Fülle von rechts auf, falls weniger als 5 Koeffizienten gegeben sind
        if len(coeffs_list) >= 1: d_ = Fraction(coeffs_list[-1])
        if len(coeffs_list) >= 2: c = Fraction(coeffs_list[-2])
        if len(coeffs_list) >= 3: b = Fraction(coeffs_list[-3])
        if len(coeffs_list) >= 4: a = Fraction(coeffs_list[-4])
        if len(coeffs_list) >= 5: e = Fraction(coeffs_list[-5])
        return cls(e,a,b,c,d_)


    def _find_roots_for_degree(self, coeffs_list_in):
        coeffs = [c for c in coeffs_list_in] # Kopie erstellen
        # Entferne führende Nullen für np.roots, aber behalte mindestens einen Koeffizienten (für konstante Fkt)
        while len(coeffs) > 1 and coeffs[0] == 0:
            coeffs.pop(0)
        
        if not coeffs: # Sollte nicht passieren, wenn wir >1 oben haben
             return ["Ungültige Koeffizienten (leer)"]

        # Spezialfälle für niedrige Grade nach Entfernung führender Nullen
        if len(coeffs) == 1: # P(x) = d (konstant)
            return ["Alle reellen Zahlen"] if coeffs[0] == 0 else []
        
        try:
            np_coeffs_float = [float(c) for c in coeffs]
            np_roots_val = np.roots(np_coeffs_float)
            real_roots_val = [r.real for r in np_roots_val if abs(r.imag) < 1e-7]
            unique_real_roots_val = []
            for r_np in sorted(real_roots_val):
                if not any(abs(r_np - ur) < 1e-7 for ur in unique_real_roots_val):
                    unique_real_roots_val.append(r_np)
            return unique_real_roots_val
        except Exception as e_np:
            print(f"Fehler bei numpy.roots mit Koeffizienten {np_coeffs_float}: {e_np}")
            return ["Numerische Lösung fehlgeschlagen"]


    def find_roots(self):
        return self._find_roots_for_degree([self.e, self.a, self.b, self.c, self.d])

    def find_critical_points_x(self):
        deriv1_c = self.derivative1_coeffs() 
        effective_deriv1_coeffs = deriv1_c[1:] # [4e, 3a, 2b, c]
        return self._find_roots_for_degree(effective_deriv1_coeffs)

    def find_inflection_points_x(self):
        deriv2_c = self.derivative2_coeffs()
        effective_deriv2_coeffs = deriv2_c[2:] # [12e, 6a, 2b]
        return self._find_roots_for_degree(effective_deriv2_coeffs)

    def expand_and_shift(self, s_shift_val):
        s = Fraction(s_shift_val)
        e0, a0, b0, c0, d0 = self.e, self.a, self.b, self.c, self.d
        e_new = e0
        a_new = a0 - 4*e0*s
        b_new = b0 - 3*a0*s + 6*e0*s**2
        c_new = c0 - 2*b0*s + 3*a0*s**2 - 4*e0*s**3
        d_new = d0 - c0*s + b0*s**2 - a0*s**3 + e0*s**4
        return Polynomial(e_new, a_new, b_new, c_new, d_new)

    def check_symmetry(self):
        # Prüft auf einfache Symmetrien
        # Achsensymmetrie zur y-Achse: P(x) = P(-x) => alle ungeraden Potenzen haben Koeffizient 0
        is_y_axis_symmetric = (self.a == 0 and self.c == 0)
        
        # Punktsymmetrie zum Ursprung: P(x) = -P(-x) => alle geraden Potenzen haben Koeffizient 0 (inkl. d)
        is_origin_symmetric = (self.e == 0 and self.b == 0 and self.d == 0)

        if is_y_axis_symmetric and is_origin_symmetric: # P(x) = 0
            return "Achsensymmetrisch zur y-Achse und Punktsymmetrisch zum Ursprung (Nullfunktion)."
        elif is_y_axis_symmetric:
            return "Achsensymmetrisch zur y-Achse."
        elif is_origin_symmetric:
            return "Punktsymmetrisch zum Ursprung."
        else:
            return "Keine einfache Symmetrie (zur y-Achse oder zum Ursprung) erkennbar."

# --- Fabrikmethoden für Basispolynome (geben Polynomial-Objekte zurück) ---
def from_saddle_point_factory(x0, D0, P_triple_prime_0, force_integer_coeffs=True):
    f_P_triple_prime_0 = Fraction(P_triple_prime_0) 
    f_x0 = Fraction(x0)
    f_D0 = Fraction(D0)
    if f_P_triple_prime_0 == 0:
        raise ValueError("P_triple_prime_0 (6a) darf nicht Null sein für kubische Funktion aus Sattelpunkt.")
    if force_integer_coeffs:
        if not (f_P_triple_prime_0.denominator == 1 and f_P_triple_prime_0.numerator % 6 == 0):
            raise ValueError("Für ganzzahlige Koeffizienten muss P_triple_prime_0 eine ganze Zahl und ein Vielfaches von 6 sein.")
    a_gen = f_P_triple_prime_0 / 6 
    a_p, b_p, c_p, d_p = a_gen, -3*f_x0*a_gen, 3*f_x0**2*a_gen, -f_x0**3*a_gen + f_D0
    if force_integer_coeffs:
        return Polynomial(0, a_p.limit_denominator(1), b_p.limit_denominator(1), c_p.limit_denominator(1), d_p.limit_denominator(1))
    else:
        return Polynomial(0, a_p, b_p, c_p, d_p)

def from_extrema_locations_factory(q1, q2, K_user, C0): 
    f_q1, f_q2, f_K_user, f_C0 = Fraction(q1), Fraction(q2), Fraction(K_user), Fraction(C0)
    if f_q1 == f_q2: raise ValueError("q1 und q2 (x-Koordinaten der Extrema) müssen verschieden sein.")
    if f_K_user == 0: raise ValueError("K_user (Skalierungsfaktor für P') darf nicht Null sein.")
    a_p = f_K_user / 3
    b_p = -f_K_user * (f_q1 + f_q2) / 2
    c_p = f_K_user * f_q1 * f_q2
    d_p = f_C0 
    return Polynomial(0, a_p, b_p, c_p, d_p) 

def from_extrema_locations_and_one_y_target_factory(q1, q2, y_at_q1, K_user, force_integer_coeffs=False): 
    f_q1, f_q2 = Fraction(q1), Fraction(q2)
    f_y_at_q1 = Fraction(y_at_q1)
    f_K_user = Fraction(K_user) 
    if f_q1 == f_q2: raise ValueError("q1 und q2 (x-Koordinaten der Extrema) müssen verschieden sein.")
    if f_K_user == 0: raise ValueError("K_user (Skalierungsfaktor für P') darf nicht Null sein.")
    term_val_at_q1 = (f_K_user/3)*f_q1**3 - (f_K_user*(f_q1+f_q2)/2)*f_q1**2 + (f_K_user*f_q1*f_q2)*f_q1
    d_p = f_y_at_q1 - term_val_at_q1
    a_p = f_K_user / 3
    b_p = -f_K_user * (f_q1 + f_q2) / 2
    c_p = f_K_user * f_q1 * f_q2
    if force_integer_coeffs:
        a_final_limited = a_p.limit_denominator(1)
        if a_p != 0 and a_final_limited == 0 : 
             raise ValueError("Mit diesen Vorgaben und 'Ganzzahlige Koeffizienten' wird der x^3-Term (a) Null. Bitte Parameter anpassen oder 'Ganzzahlige Koeffizienten' deaktivieren.")
        return Polynomial(0, a_final_limited,b_p.limit_denominator(1),c_p.limit_denominator(1),d_p.limit_denominator(1))
    return Polynomial(0, a_p, b_p, c_p, d_p)

def from_extrema_locations_and_values_factory(q1, q2, y1, y2, force_integer_coeffs=False): 
    f_q1, f_q2 = Fraction(q1), Fraction(q2)
    f_y1, f_y2 = Fraction(y1), Fraction(y2)
    if f_q1 == f_q2:
        raise ValueError("q1 und q2 (x-Koordinaten der Extrema) müssen verschieden sein.")
    def term_T_val(x_val, f_q1_in, f_q2_in):
        x = Fraction(x_val)
        return (x**3 / 3) - ((f_q1_in + f_q2_in) / 2) * x**2 + (f_q1_in * f_q2_in) * x
    T_at_q1 = term_T_val(f_q1, f_q1, f_q2)
    T_at_q2 = term_T_val(f_q2, f_q1, f_q2)
    if T_at_q2 == T_at_q1:
        if f_y1 == f_y2: raise ValueError("K_user unbestimmt: T(q1)=T(q2) und y1=y2. Dies führt zu einem konstanten Polynom oder unendlich vielen Lösungen.")
        else: raise ValueError("Keine Lösung für K_user: T(q1)=T(q2) aber y1!=y2.")
    K_user = (f_y2 - f_y1) / (T_at_q2 - T_at_q1)
    
    if K_user == 0: 
        raise ValueError("Die angegebenen y-Werte für die Extrema führen zu K_user = 0 (und damit a=0). Ein kubisches Polynom mit diesen Eigenschaften ist nicht eindeutig oder nicht existent.")

    D_const = f_y1 - K_user * T_at_q1
    a_p = K_user / 3
    b_p = -K_user * (f_q1 + f_q2) / 2
    c_p = K_user * f_q1 * f_q2
    d_p = D_const
    if force_integer_coeffs:
        a_final_limited = a_p.limit_denominator(1)
        if a_p != 0 and a_final_limited == 0: 
            raise ValueError("Mit diesen Extrema-Vorgaben und 'Ganzzahlige Koeffizienten' wird der x^3-Term (a) Null. Bitte Parameter anpassen oder 'Ganzzahlige Koeffizienten' deaktivieren.")
        return Polynomial(0, a_final_limited,b_p.limit_denominator(1),c_p.limit_denominator(1),d_p.limit_denominator(1))
    else:
        return Polynomial(0, a_p, b_p, c_p, d_p)

def from_inflection_point_and_delta_q_factory(xw, delta_q, K_user, C0): 
    f_xw,f_delta_q,f_K_user,f_C0 = Fraction(xw),Fraction(delta_q),Fraction(K_user),Fraction(C0)
    if f_delta_q == 0: raise ValueError("delta_q (Abstand zu Extrema) muss ungleich Null sein.")
    if f_K_user == 0: raise ValueError("K_user (Skalierungsfaktor für P') darf nicht Null sein.")
    a_p = f_K_user / 3
    b_p = -f_K_user * f_xw
    c_p = f_K_user * (f_xw**2 - f_delta_q**2)
    d_p = f_C0
    return Polynomial(0, a_p, b_p, c_p, d_p)

def from_integer_inflection_real_extrema_factory(xw, D0, a_coeff, C1_coeff): 
    f_xw,f_D0,f_a_coeff,f_C1_coeff = Fraction(xw),Fraction(D0),Fraction(a_coeff),Fraction(C1_coeff)
    if f_a_coeff == 0: raise ValueError("a_coeff (Koeffizient von x^3) darf nicht Null sein.")
    if f_C1_coeff != 0 and (f_C1_coeff * f_a_coeff) >= 0: 
        raise ValueError("Für reelle Extrema müssen a_coeff und C1_coeff unterschiedliche Vorzeichen haben.")
    A_ = f_a_coeff 
    B_ = -3*f_a_coeff*f_xw
    C_ = 3*f_a_coeff*f_xw**2 + f_C1_coeff
    D_ = -f_a_coeff*f_xw**3 - f_C1_coeff*f_xw + f_D0
    return Polynomial(0, A_,B_,C_,D_)

def from_zero_and_two_integer_roots_factory(r2, r3, scale_factor_S): 
    f_r2,f_r3,f_S = Fraction(r2),Fraction(r3),Fraction(scale_factor_S)
    if f_S == 0: raise ValueError("scale_factor_S darf nicht Null sein.")
    a_p = f_S
    b_p = -f_S*(f_r2+f_r3)
    c_p = f_S*f_r2*f_r3
    d_p = Fraction(0)
    return Polynomial(0, a_p,b_p,c_p,d_p)
    
def from_quadratic_vertex_and_scale_factory(h, k, scale_factor_sq, force_integer_coeffs=False): 
    f_h, f_k, f_scale = Fraction(h), Fraction(k), Fraction(scale_factor_sq)
    if f_scale == 0: raise ValueError("scale_factor_sq (Koeffizient b für x^2) darf nicht Null sein.")
    b_p = f_scale
    c_p = -2*f_scale*f_h
    d_p = f_scale*f_h**2 + f_k
    if force_integer_coeffs:
        b_final_limited = b_p.limit_denominator(1)
        if b_p != 0 and b_final_limited == 0: 
            raise ValueError("Mit diesen Scheitelpunkt-Vorgaben und 'Ganzzahlige Koeffizienten' wird der x^2-Term (b) Null. Bitte Parameter anpassen oder 'Ganzzahlige Koeffizienten' deaktivieren.")
        return Polynomial(0, 0, b_final_limited, c_p.limit_denominator(1), d_p.limit_denominator(1))
    else:
        return Polynomial(0, 0, b_p, c_p, d_p)

def from_linear_root_and_y_intercept_factory(root_x, y_intercept): 
    f_root_x, f_y_intercept = Fraction(root_x), Fraction(y_intercept)
    d_p = f_y_intercept
    c_p = Fraction(0)
    if f_root_x == 0:
        if f_y_intercept != 0: 
            raise ValueError("Wenn die Nullstelle x=0 ist, muss der y-Achsenabschnitt auch 0 sein (P(0)=0).")
        else: 
              c_p = Fraction(0) 
    else: 
        c_p = -f_y_intercept / f_root_x
    return Polynomial(0, 0, 0, c_p, d_p)

# --- Fabrikmethoden für Biquadratische Basispolynome ---
def from_outer_extrema_and_central_extremum_biquad_factory(x_outer_ext, y_outer_ext, y_central_ext, force_integer_coeffs=False):
    f_x_e = Fraction(x_outer_ext)
    f_y_e = Fraction(y_outer_ext)
    f_y_0 = Fraction(y_central_ext)
    if f_x_e == 0:
        raise ValueError("x-Koordinate des äußeren Extremums (x_e) darf nicht Null sein.")
    D_const = f_y_0
    if f_x_e**4 == 0 : 
         raise ValueError("x_e^4 darf nicht Null sein (x_e darf nicht Null sein).")
    A_p = (D_const - f_y_e) / (f_x_e**4)
    B_p = -2 * A_p * (f_x_e**2)
    if force_integer_coeffs:
        A_f_limited = A_p.limit_denominator(1)
        if A_p != 0 and A_f_limited == 0: 
            raise ValueError("Koeffizient A wird bei Ganzzahligkeit zu Null. Bitte Parameter anpassen oder 'Ganzzahlige Koeffizienten' deaktivieren.")
        return Polynomial(A_f_limited, 0, B_p.limit_denominator(1), 0, D_const.limit_denominator(1)) # e,a,b,c,d
    return Polynomial(A_p, 0, B_p, 0, D_const)

def from_inflection_points_and_central_extremum_biquad_factory(x_wp, y_wp, y_central_ext, force_integer_coeffs=False):
    f_x_wp = Fraction(x_wp)
    f_y_wp = Fraction(y_wp)
    f_y_0 = Fraction(y_central_ext)
    if f_x_wp == 0:
        raise ValueError("x-Koordinate des Wendepunkts (x_wp) darf nicht Null sein für diese Konstruktion.")
    D_const = f_y_0 
    denom_A = 5 * f_x_wp**4
    if denom_A == 0: 
        raise ValueError("5 * x_wp^4 darf nicht Null sein (x_wp darf nicht Null sein).")
    A_p = (D_const - f_y_wp) / denom_A
    B_p = -6 * A_p * (f_x_wp**2)
    if force_integer_coeffs:
        A_f_limited = A_p.limit_denominator(1)
        if A_p != 0 and A_f_limited == 0:
            raise ValueError("Koeffizient A wird bei Ganzzahligkeit zu Null. Bitte Parameter anpassen oder 'Ganzzahlige Koeffizienten' deaktivieren.")
        return Polynomial(A_f_limited, 0, B_p.limit_denominator(1), 0, D_const.limit_denominator(1))
    return Polynomial(A_p, 0, B_p, 0, D_const)

def direct_coeffs_biquad_factory(A, B, D_val):
    return Polynomial(A,0,B,0,D_val)


class ExportDialog(tk.Toplevel):
    def __init__(self, parent, title, text_content):
        super().__init__(parent)
        self.title(title)
        lines = text_content.count('\n') + 1
        width_char = 80 
        height_char = min(max(lines + 5, 15), 35) 
        
        width_px = width_char * 8 
        height_px = height_char * 18
        self.geometry(f"{width_px}x{height_px}")


        main_frame = ttk.Frame(self, padding=10)
        main_frame.pack(expand=True, fill=tk.BOTH)

        text_area = tk.Text(main_frame, wrap=tk.WORD, font=("Consolas", 10), relief=tk.SOLID, borderwidth=1)
        text_area.insert(tk.END, text_content)
        text_area.config(state=tk.DISABLED)
        
        v_scrollbar = ttk.Scrollbar(main_frame, orient=tk.VERTICAL, command=text_area.yview)
        text_area['yscrollcommand'] = v_scrollbar.set
        
        v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        text_area.pack(side=tk.LEFT, expand=True, fill=tk.BOTH)

        button_frame = ttk.Frame(self)
        button_frame.pack(pady=5)

        copy_button = ttk.Button(button_frame, text="In Zwischenablage kopieren", command=lambda: self.copy_to_clipboard(text_content))
        copy_button.pack(side=tk.LEFT, padx=5)
        
        close_button = ttk.Button(button_frame, text="Schließen", command=self.destroy)
        close_button.pack(side=tk.LEFT, padx=5)

        self.protocol("WM_DELETE_WINDOW", self.destroy)
        self.transient(parent) 
        self.grab_set() 
        self.focus_set()
        self.wait_window()

    def copy_to_clipboard(self, text):
        self.clipboard_clear()
        self.clipboard_append(text)
        messagebox.showinfo("Kopiert", "Text wurde in die Zwischenablage kopiert.", parent=self)


class PolynomialFeatures:
    def __init__(self, poly_obj):
        self.poly = poly_obj
        self.precision = 3 

    def _format_number_list(self, data_list):
        if data_list is None or (isinstance(data_list, list) and not data_list):
            return "Keine"
        if isinstance(data_list, list):
            if data_list and isinstance(data_list[0], str): 
                return data_list[0]
            try:
                return ", ".join([f"{float(x):.{self.precision}f}" for x in data_list])
            except (ValueError, TypeError):
                return "Formatierungsfehler"
        try: 
            return f"{float(data_list):.{self.precision}f}"
        except (ValueError, TypeError):
            return "Formatierungsfehler"

    def get_polynomial_string(self):
        return f"P(x) = {str(self.poly)}"

    def get_derivatives_text(self):
        p_prime_obj = Polynomial.from_coeffs(self.poly.derivative1_coeffs()[1:])
        p_double_prime_obj = Polynomial.from_coeffs(self.poly.derivative2_coeffs()[2:])
        p_triple_prime_coeffs = p_double_prime_obj.derivative1_coeffs() 
        p_triple_prime_obj = Polynomial.from_coeffs(p_triple_prime_coeffs[1:])


        text = "1. Ableitungen:\n"
        text += f"   P'(x)  = {str(p_prime_obj)}\n"
        text += f"   P''(x) = {str(p_double_prime_obj)}\n"
        text += f"   P'''(x) = {str(p_triple_prime_obj)}\n"
        return text

    def get_roots_text(self):
        roots = self.poly.find_roots()
        return f"2. Nullstellen:\n   x_N = {self._format_number_list(roots)}\n"

    def get_extrema_text(self):
        crit_x = self.poly.find_critical_points_x()
        text = "3. Extrema (Kritische Punkte):\n"
        if not crit_x or (isinstance(crit_x, list) and crit_x and isinstance(crit_x[0], str)):
            text += f"   x_E = {self._format_number_list(crit_x)}\n"
            text += "   Keine diskreten Extrema gefunden.\n"
            return text

        text += f"   x_E = {self._format_number_list(crit_x)}\n"
        
        p_double_prime_obj = Polynomial.from_coeffs(self.poly.derivative2_coeffs()[2:])
        p_triple_prime_obj = Polynomial.from_coeffs(p_double_prime_obj.derivative1_coeffs()[1:])


        extrema_details = []
        for x_val_frac in crit_x: # x_val ist hier Fraction oder float
            x_val = float(x_val_frac) # Für Auswertung in Ableitungen
            y_val = self.poly.evaluate(x_val_frac) # Für y-Wert exakt
            p_double_val = p_double_prime_obj.evaluate(x_val) # Numerisch für Klassifizierung
            
            point_str = f"({self._format_number_list(x_val)} | {self._format_number_list(y_val)})"
            
            if abs(p_double_val) < 1e-7: 
                p_triple_val = p_triple_prime_obj.evaluate(x_val)
                if abs(p_triple_val) < 1e-7: 
                    extrema_details.append(f"   - {point_str}: P''(x_E)=0, P'''(x_E)=0 (Weitere Untersuchung nötig, evtl. Flachpunkt)")
                else: 
                    extrema_details.append(f"   - Sattelpunkt (Terrassenpunkt) bei {point_str} (P''(x_E)=0, P'''(x_E)≠0)")
            elif p_double_val > 0:
                extrema_details.append(f"   - Tiefpunkt (Minimum) bei {point_str} (P''(x_E) > 0)")
            else: 
                extrema_details.append(f"   - Hochpunkt (Maximum) bei {point_str} (P''(x_E) < 0)")
        
        if extrema_details:
            text += "\n".join(extrema_details) + "\n"
        else:
            text += "   Keine klassifizierbaren Extrema gefunden.\n"
            
        return text

    def get_inflection_points_text(self):
        infl_x = self.poly.find_inflection_points_x()
        text = "4. Wendepunkte:\n"
        if not infl_x or (isinstance(infl_x, list) and infl_x and isinstance(infl_x[0], str)):
            text += f"   x_W = {self._format_number_list(infl_x)}\n"
            return text

        text += f"   x_W = {self._format_number_list(infl_x)}\n"
        
        wp_details = []
        for x_val_frac in infl_x:
            x_val = float(x_val_frac)
            y_val = self.poly.evaluate(x_val_frac)
            wp_details.append(f"   - Wendepunkt bei ({self._format_number_list(x_val)} | {self._format_number_list(y_val)})")
        
        if wp_details:
            text += "\n".join(wp_details) + "\n"
        else:
            text += "   Keine Wendepunkte gefunden.\n"
        return text

    def get_symmetry_text(self):
        return f"5. Symmetrie:\n   {self.poly.check_symmetry()}\n"

    def get_monotonicity_text(self):
        text = "6. Monotonieverhalten:\n"
        crit_points_raw = self.poly.find_critical_points_x()
        
        # Filtere nicht-numerische Werte aus kritischen Punkten (z.B. "Alle reellen Zahlen")
        crit_points = [cp for cp in crit_points_raw if isinstance(cp, (int, float, Fraction))]

        if not crit_points:
            p_prime_obj = Polynomial.from_coeffs(self.poly.derivative1_coeffs()[1:])
            # Test an einem beliebigen Punkt, wenn keine kritischen Punkte vorhanden sind
            # Dies ist relevant, wenn P'(x) keine Nullstellen hat (z.B. P'(x) = x^2 + 1)
            # oder wenn P'(x) konstant ist.
            is_const_deriv = (p_prime_obj.e == 0 and p_prime_obj.a == 0 and p_prime_obj.b == 0 and p_prime_obj.c ==0)

            if is_const_deriv: # P'(x) = d'
                if p_prime_obj.d > 0: text += "   Streng monoton steigend in R.\n"
                elif p_prime_obj.d < 0: text += "   Streng monoton fallend in R.\n"
                else: text += "   Konstant in R (Steigung 0).\n" # P'(x) = 0
            else: # P'(x) hat keine Nullstellen, aber ist nicht konstant 0
                test_val = p_prime_obj.evaluate(0) 
                if test_val > 1e-7: text += "   Streng monoton steigend in R.\n"
                elif test_val < -1e-7: text += "   Streng monoton fallend in R.\n"
                else: text += "   Monotonieverhalten komplex (P'(x) könnte überall 0 sein oder oszillieren ohne Nullstellen).\n"
            return text

        sorted_crit_points = sorted(list(set(float(cp) for cp in crit_points)))
        intervals = []
        if not sorted_crit_points: # Sollte nach obiger Prüfung nicht mehr eintreten, aber als Fallback
            text += "   Keine kritischen Punkte für Intervallanalyse gefunden.\n"
            return text

        if len(sorted_crit_points) == 1:
            intervals.append((float('-inf'), sorted_crit_points[0]))
            intervals.append((sorted_crit_points[0], float('inf')))
        else:
            intervals.append((float('-inf'), sorted_crit_points[0]))
            for i in range(len(sorted_crit_points) - 1):
                intervals.append((sorted_crit_points[i], sorted_crit_points[i+1]))
            intervals.append((sorted_crit_points[-1], float('inf')))

        p_prime_obj = Polynomial.from_coeffs(self.poly.derivative1_coeffs()[1:])
        for interval in intervals:
            lower, upper = interval
            if lower == float('-inf') and upper == float('inf'): test_x = 0
            elif lower == float('-inf'): test_x = upper - 1
            elif upper == float('inf'): test_x = lower + 1
            else: test_x = (lower + upper) / 2
            
            epsilon = 1e-5 # Kleine Verschiebung, falls Testpunkt exakt auf kritischem Punkt landet
            # Prüfe, ob test_x zu nah an einem kritischen Punkt ist
            is_too_close = False
            for cp in sorted_crit_points:
                if abs(test_x - cp) < epsilon:
                    is_too_close = True
                    break
            if is_too_close: # Versuche, den Testpunkt anzupassen
                if lower != float('-inf') and upper != float('inf'):
                    test_x = lower + (upper - lower) * 0.25 # Ein anderer Punkt im Intervall
                    if any(abs(test_x - cp) < epsilon for cp in sorted_crit_points):
                         test_x = lower + (upper - lower) * 0.75
                elif lower == float('-inf'):
                     test_x = upper - abs(upper*0.1 if upper !=0 else 1) - epsilon # Weiter weg von der oberen Grenze
                elif upper == float('inf'):
                     test_x = lower + abs(lower*0.1 if lower !=0 else 1) + epsilon # Weiter weg von der unteren Grenze


            p_prime_at_test_x = p_prime_obj.evaluate(test_x)
            interval_str = f"({self._format_number_list(lower) if lower != float('-inf') else '-∞'}; {self._format_number_list(upper) if upper != float('inf') else '∞'})"
            if p_prime_at_test_x > 1e-7: 
                text += f"   Streng monoton steigend in {interval_str}\n"
            elif p_prime_at_test_x < -1e-7: 
                text += f"   Streng monoton fallend in {interval_str}\n"
            else: 
                text += f"   P'(x) ≈ 0 in {interval_str} (evtl. konstant oder Sattelpunktbereich)\n"
        return text

    def get_curvature_text(self):
        text = "7. Krümmungsverhalten:\n"
        infl_points_raw = self.poly.find_inflection_points_x()
        infl_points = [ip for ip in infl_points_raw if isinstance(ip, (int, float, Fraction))]


        if not infl_points:
            p_double_prime_obj = Polynomial.from_coeffs(self.poly.derivative2_coeffs()[2:])
            is_const_double_deriv = (p_double_prime_obj.e == 0 and p_double_prime_obj.a == 0 and p_double_prime_obj.b == 0 and p_double_prime_obj.c ==0)
            if is_const_double_deriv:
                if p_double_prime_obj.d > 0: text += "   Konvex (linksgekrümmt) in R.\n"
                elif p_double_prime_obj.d < 0: text += "   Konkav (rechtsgekrümmt) in R.\n"
                else: text += "   P''(x) = 0 in R (keine Krümmung, z.B. lineare Funktion).\n"
            else: # P''(x) hat keine Nullstellen, aber ist nicht konstant 0
                test_val = p_double_prime_obj.evaluate(0)
                if test_val > 1e-7: text += "   Konvex (linksgekrümmt) in R.\n"
                elif test_val < -1e-7: text += "   Konkav (rechtsgekrümmt) in R.\n"
                else: text += "   Krümmungsverhalten komplex (P''(x) könnte überall 0 sein oder oszillieren ohne Nullstellen).\n"
            return text

        sorted_infl_points = sorted(list(set(float(ip) for ip in infl_points)))
        intervals = []
        if not sorted_infl_points:
            text += "   Keine Wendepunkte für Intervallanalyse gefunden.\n"
            return text

        if len(sorted_infl_points) == 1:
            intervals.append((float('-inf'), sorted_infl_points[0]))
            intervals.append((sorted_infl_points[0], float('inf')))
        else:
            intervals.append((float('-inf'), sorted_infl_points[0]))
            for i in range(len(sorted_infl_points) - 1):
                intervals.append((sorted_infl_points[i], sorted_infl_points[i+1]))
            intervals.append((sorted_infl_points[-1], float('inf')))

        p_double_prime_obj = Polynomial.from_coeffs(self.poly.derivative2_coeffs()[2:])
        for interval in intervals:
            lower, upper = interval
            if lower == float('-inf') and upper == float('inf'): test_x = 0
            elif lower == float('-inf'): test_x = upper - 1
            elif upper == float('inf'): test_x = lower + 1
            else: test_x = (lower + upper) / 2
            
            epsilon = 1e-5
            while any(abs(test_x - ip) < epsilon for ip in sorted_infl_points):
                test_x += epsilon 

            p_double_at_test_x = p_double_prime_obj.evaluate(test_x)
            interval_str = f"({self._format_number_list(lower) if lower != float('-inf') else '-∞'}; {self._format_number_list(upper) if upper != float('inf') else '∞'})"
            if p_double_at_test_x > 1e-7:
                text += f"   Konvex (linksgekrümmt) in {interval_str}\n"
            elif p_double_at_test_x < -1e-7:
                text += f"   Konkav (rechtsgekrümmt) in {interval_str}\n"
            else: 
                text += f"   P''(x) ≈ 0 in {interval_str} (genauere Untersuchung nötig)\n"
        return text


    def get_full_discussion(self):
        discussion = "Kurvendiskussion für:\n"
        discussion += self.get_polynomial_string() + "\n\n"
        # Definitionsbereich wird weggelassen
        discussion += self.get_symmetry_text() + "\n"
        discussion += self.get_derivatives_text() + "\n"
        discussion += self.get_roots_text() + "\n"
        discussion += self.get_extrema_text() + "\n"
        discussion += self.get_inflection_points_text() + "\n"
        discussion += self.get_monotonicity_text() + "\n"
        discussion += self.get_curvature_text()
        return discussion


class PolynomialAppGUI:
    MAX_PARAM_ROWS = 6 

    def __init__(self, master, initial_width=1440, initial_height=970): # Standardgröße angepasst
        self.master = master
        master.title("Polynomfunktionen Generator")
        master.geometry(f"{initial_width}x{initial_height}") 

        self.style = ttk.Style()
        self.style.theme_use('clam') 

        self.paned_window = ttk.PanedWindow(master, orient=tk.HORIZONTAL)
        self.paned_window.pack(expand=True, fill=tk.BOTH, padx=10, pady=10)

        controls_outer_frame = ttk.Frame(self.paned_window, padding="5")
        self.paned_window.add(controls_outer_frame, weight=1) 

        plot_outer_frame = ttk.LabelFrame(self.paned_window, text="Graphische Darstellung", padding="5")
        self.paned_window.add(plot_outer_frame, weight=2) 

        # Frame für String-Eingabe
        string_input_outer_frame = ttk.Frame(controls_outer_frame)
        string_input_outer_frame.pack(pady=(5,0), fill=tk.X) 

        string_input_frame = ttk.LabelFrame(string_input_outer_frame, text="Polynom direkt eingeben", padding="10")
        string_input_frame.pack(pady=0, fill=tk.X) 
        
        self.poly_string_input_var = tk.StringVar()
        ttk.Label(string_input_frame, text="P(x) =").pack(side=tk.LEFT, padx=(0,5))
        self.poly_string_entry = ttk.Entry(string_input_frame, textvariable=self.poly_string_input_var, width=45) 
        self.poly_string_entry.pack(side=tk.LEFT, padx=0, expand=True, fill=tk.X)
        
        self.clear_string_button = ttk.Button(string_input_frame, text="X", command=lambda: self.poly_string_input_var.set(""), width=3)
        self.clear_string_button.pack(side=tk.LEFT, padx=(2,5))

        self.parse_string_button = ttk.Button(string_input_frame, text="Übernehmen & Analysieren", command=self.apply_parsed_polynomial)
        self.parse_string_button.pack(side=tk.LEFT, padx=(5,0))


        poly_type_frame = ttk.LabelFrame(controls_outer_frame, text="Polynomtyp", padding="10")
        poly_type_frame.pack(pady=5, fill=tk.X)
        self.poly_type_var = tk.StringVar(value="Kubisch")
        poly_types = ["Kubisch", "Biquadratisch", "Quadratisch", "Linear"] 
        ttk.Label(poly_type_frame, text="Typ wählen:").pack(side=tk.LEFT, padx=5)
        self.poly_type_menu = ttk.Combobox(poly_type_frame, textvariable=self.poly_type_var,
                                           values=poly_types, state="readonly", width=20)
        self.poly_type_menu.pack(side=tk.LEFT, padx=5)
        self.poly_type_menu.bind("<<ComboboxSelected>>", self.update_ui_for_poly_type)

        self.construction_frame = ttk.LabelFrame(controls_outer_frame, text="Primäre Konstruktionseigenschaft", padding="10")
        
        self.construction_method_var = tk.StringVar()
        self.cubic_construction_options = {
            "Aus Koeffizienten (a,b,c,d)": "direct_coeffs_cubic", 
            "Sattelpunkt definieren": "from_saddle_point",
            "Extrema-Lagen (x-Werte) definieren": "from_extrema_locations",
            "Extrema (x-Werte) & y-Ziel für q1": "from_extrema_locations_and_one_y_target", 
            "Extrema (x- und y-Werte) definieren": "from_extrema_locations_and_values", 
            "WP & Abstand zu Extrema definieren": "from_inflection_point_and_delta_q",
            "WP (x-Wert), Extrema reell": "from_integer_inflection_real_extrema",
            "Nullstellen (0, r2, r3) definieren": "from_zero_and_two_integer_roots"
        }
        self.biquadratic_construction_options = { 
            "Aus Koeffizienten (A,B,D) für Ax^4+Bx^2+D": "direct_coeffs_biquad", 
            "Äußere Extrema (x_e, y_e) & Zentrales Extremum (y_0)": "from_outer_extrema_and_central_extremum_biquad",
            "Wendepunkte (x_wp, y_wp) & Zentrales Extremum (y_0)": "from_inflection_points_and_central_extremum_biquad"
        }
        self.quadratic_construction_options = {
            "Aus Koeffizienten (b,c,d) für bx^2+cx+d": "direct_coeffs_quad", 
            "Aus Scheitelpunkt und Skalierungsfaktor": "from_quadratic_vertex_and_scale"
        }
        self.linear_construction_options = {
            "Aus Koeffizienten (c,d) für cx+d": "direct_coeffs_lin", 
            "Aus Nullstelle und y-Achsenabschnitt": "from_linear_root_and_y_intercept"
        }
        
        self.construction_method_label = ttk.Label(self.construction_frame, text="Konstruktion durch:")
        self.construction_method_menu = ttk.Combobox(self.construction_frame, 
                                                     textvariable=self.construction_method_var,
                                                     state="readonly", width=40) 
        self.construction_method_menu.bind("<<ComboboxSelected>>", self._update_construction_params_display) 
        
        self.modifier_frame = ttk.LabelFrame(controls_outer_frame, text="Globale Modifikatoren", padding="10") 
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
        
        self.input_params_frame = ttk.LabelFrame(controls_outer_frame, text="Parameter", padding="10")
        self.input_params_frame.pack(pady=5, fill=tk.X, expand=False) 
        self.input_params_frame.columnconfigure(0, weight=0, minsize=200)  
        self.input_params_frame.columnconfigure(1, weight=0, minsize=150)  
        self.input_params_frame.columnconfigure(2, weight=1, minsize=100)  
        
        self.param_row_widgets = []
        self.param_row_min_height = 30 
        for i in range(self.MAX_PARAM_ROWS):
            self.input_params_frame.rowconfigure(i, minsize=self.param_row_min_height)
            lbl = ttk.Label(self.input_params_frame, text="")
            entry_var = tk.StringVar()
            entry = ttk.Entry(self.input_params_frame, textvariable=entry_var, width=20) 
            hint_lbl = ttk.Label(self.input_params_frame, text="", foreground="blue", font=("TkDefaultFont", 8))
            self.param_row_widgets.append({'label': lbl, 'entry_var': entry_var, 'entry': entry, 'hint': hint_lbl})
            
        self.param_entries = {} 
        self.analyzed_poly = None 

        button_frame = ttk.Frame(controls_outer_frame) 
        button_frame.pack(pady=10)

        self.generate_button = ttk.Button(button_frame, text="Polynom generieren und analysieren", command=self.generate_and_display) # self.generate_button
        self.generate_button.pack(side=tk.LEFT, padx=5)
        
        self.flip_extrema_button = ttk.Button(button_frame, text="Extrema-Typ umkehren (K)", command=self.flip_extrema_type)
        
        lucky_button = ttk.Button(button_frame, text="Ich fühle mich glücklich!", command=self.i_feel_lucky)
        lucky_button.pack(side=tk.LEFT, padx=5)

        output_frame = ttk.LabelFrame(controls_outer_frame, text="Ergebnisse", padding="10") 
        output_frame.pack(pady=5, expand=True, fill=tk.BOTH)
        output_frame.columnconfigure(1, weight=1) 

        self.output_fields = {}
        output_labels_texts = { 
            "P(x)": "Polynom P(x):", "P'(x)": "Ableitung P'(x):", "P''(x)": "Ableitung P''(x):",
            "roots": "Nullstellen P(x):",
            "crit_x": "Kritische Stellen x:", "crit_y": "Y-Werte an krit. Stellen:",
            "infl_x": "Wendestelle(n) x:", "infl_y": "Y-Werte an Wendepunkt(en):",
            "saddle_info": "Sattelpunkt Info:"
        }
        row_idx = 0
        for key, text in output_labels_texts.items():
            ttk.Label(output_frame, text=text).grid(row=row_idx, column=0, sticky=tk.NW, padx=5, pady=2)
            text_widget = tk.Text(output_frame, height=1, width=45, wrap=tk.WORD, relief=tk.SUNKEN, borderwidth=1, font=("Consolas", 10)) 
            if key in ["P(x)", "P'(x)", "P''(x)"]: text_widget.config(height=2)
            text_widget.grid(row=row_idx, column=1, sticky=tk.EW, padx=5, pady=2)
            text_widget.config(state=tk.DISABLED) 
            self.output_fields[key] = text_widget
            row_idx += 1
        
        self.export_discussion_button = ttk.Button(output_frame, text="Kurvendiskussion exportieren", command=self.export_curve_discussion, state=tk.DISABLED)
        self.export_discussion_button.grid(row=row_idx, column=0, columnspan=2, pady=10)

        # Statusleiste HINZUGEFÜGT
        self.status_label_var = tk.StringVar()
        status_bar = ttk.Label(controls_outer_frame, textvariable=self.status_label_var, relief=tk.SUNKEN, anchor=tk.W, padding=2)
        status_bar.pack(side=tk.BOTTOM, fill=tk.X)
        self.update_status("Bereit.")


        self.fig = Figure(figsize=(6, 4.5), dpi=100) 
        self.plot_ax = self.fig.add_subplot(111)
        self.plot_ax.set_xlabel("x")
        self.plot_ax.set_ylabel("P(x)") 
        self.plot_ax.grid(True, linestyle=':', alpha=0.7)
        self.plot_ax.axhline(0, color='black', linewidth=0.5)
        self.plot_ax.axvline(0, color='black', linewidth=0.5)

        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_outer_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        
        self.toolbar = NavigationToolbar2Tk(self.canvas, plot_outer_frame)
        self.toolbar.update()
        self.toolbar.pack(side=tk.BOTTOM, fill=tk.X) 
        self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.canvas.draw() 
        self.update_ui_for_poly_type() 

    def update_status(self, message):
        self.status_label_var.set(message)

    def on_modifier_toggle(self):
        self._apply_param_hints()

    def update_ui_for_poly_type(self, event=None):
        poly_type = self.poly_type_var.get()
        self.clear_param_input_fields_display() 
        self.flip_extrema_button.pack_forget() 

        current_options = {}
        frame_text_prefix = "Primäre Konstruktionseigenschaft"
        
        self.x_shift_specific_param_name = "x_shift_val" 

        if poly_type == "Kubisch":
            current_options = self.cubic_construction_options
            frame_text_prefix = "Konstruktionseigenschaft (Kubisch)"
            self.modifier_frame.pack(pady=5, fill=tk.X, before=self.input_params_frame)
        elif poly_type == "Biquadratisch": 
            current_options = self.biquadratic_construction_options
            frame_text_prefix = "Konstruktionsmethode (Biquadratisch)"
            self.modifier_frame.pack(pady=5, fill=tk.X, before=self.input_params_frame)
        elif poly_type == "Quadratisch":
            current_options = self.quadratic_construction_options
            frame_text_prefix = "Konstruktionsmethode (Quadratisch)"
            self.modifier_frame.pack(pady=5, fill=tk.X, before=self.input_params_frame)
        elif poly_type == "Linear":
            current_options = self.linear_construction_options
            frame_text_prefix = "Konstruktionsmethode (Linear)"
            self.modifier_frame.pack_forget() 

        self.construction_frame.config(text=frame_text_prefix)
        self.construction_method_label.pack(side=tk.LEFT, padx=5)
        self.construction_method_menu.config(values=list(current_options.keys()))
        self.construction_method_menu.pack(side=tk.LEFT, padx=5, expand=True, fill=tk.X)
        self.construction_frame.pack(pady=5, fill=tk.X, before=self.input_params_frame) 
        
        current_selection = self.construction_method_var.get()
        if not current_selection or current_selection not in current_options:
            if list(current_options.keys()): 
                 self.construction_method_menu.current(0) 
        
        self._update_construction_params_display()


    def clear_param_input_fields_display(self):
        for i in range(self.MAX_PARAM_ROWS):
            row = self.param_row_widgets[i]
            row['label'].grid_forget()
            row['entry'].grid_forget()
            row['hint'].grid_forget()
            row['label'].config(text="")
            row['entry_var'].set("") 
            row['hint'].config(text="")
        self.param_entries.clear() 

    def _configure_param_entry_row(self, row_index, param_name, label_text, default_value=""):
        if row_index >= self.MAX_PARAM_ROWS:
            print(f"Warnung: Zeile {row_index} überschreitet MAX_PARAM_ROWS.")
            return row_index 
        row = self.param_row_widgets[row_index]
        row['label'].config(text=label_text)
        row['entry_var'].set(str(default_value)) 
        row['label'].grid(row=row_index, column=0, sticky=tk.W, padx=5, pady=3)
        row['entry'].grid(row=row_index, column=1, sticky=tk.W, padx=5, pady=3) 
        row['hint'].grid(row=row_index, column=2, sticky=tk.W, padx=5, pady=3) 
        self.param_entries[param_name] = row['entry_var'] 
        return row_index + 1
    
    def _add_x_shift_parameter_if_applicable(self, current_row):
        poly_type = self.poly_type_var.get()
        if poly_type in ["Kubisch", "Biquadratisch"]:
            current_row = self._configure_param_entry_row(current_row, self.x_shift_specific_param_name, "X-Verschiebung (s):", "0")
        return current_row


    def _update_construction_params_display(self, event=None): 
        poly_type = self.poly_type_var.get()
        if poly_type == "Kubisch":
            self._update_cubic_construction_params_create_fields()
        elif poly_type == "Biquadratisch": 
            self._update_biquadratic_construction_params_create_fields()
        elif poly_type == "Quadratisch":
            self._update_quadratic_construction_params_create_fields()
        elif poly_type == "Linear":
            self._update_linear_construction_params_create_fields()
        
    def _apply_param_hints(self):
        poly_type = self.poly_type_var.get()
        construction_key = self.construction_method_var.get()
        force_int_active = self.force_int_coeffs_var.get()
        aim_int_y_active = self.aim_int_y_var.get()

        for i in range(self.MAX_PARAM_ROWS):
            if self.param_row_widgets[i]['label'].winfo_ismapped(): 
                 self.param_row_widgets[i]['hint'].config(text="")

        if poly_type == "Kubisch":
            method_name = self.cubic_construction_options.get(construction_key)
            if method_name == "from_saddle_point":
                if "P_triple_prime_0" in self.param_entries: self.param_row_widgets[2]['hint'].config(text="Ganze Zahl & Vielfaches von 6" if force_int_active else "Kann Bruch sein (z.B. 1/6)")
                if aim_int_y_active and "D0" in self.param_entries: self.param_row_widgets[1]['hint'].config(text="Ganzzahlig für ganzz. y am SP")
            elif method_name == "from_extrema_locations":
                if aim_int_y_active and "C0" in self.param_entries: self.param_row_widgets[3]['hint'].config(text="Ganzzahlig für ganzz. y-Extrema (beeinflusst d)")
                if "K_user" in self.param_entries: self.param_row_widgets[2]['hint'].config(text="K für P'. Vorzeichen bestimmt Min/Max-Reihenfolge. Vielfaches von 6 für ganzz. Koeff.")
            elif method_name == "from_extrema_locations_and_one_y_target":
                if aim_int_y_active and "y_at_q1_onetarget" in self.param_entries: self.param_row_widgets[2]['hint'].config(text="Ganzzahliger Ziel-y-Wert")
                if "K_user_onetarget" in self.param_entries: self.param_row_widgets[3]['hint'].config(text="K für P'. Vielfaches von 6 für ganzz. Koeff.")
                if force_int_active and "q1_onetarget" in self.param_entries : self.param_row_widgets[0]['hint'].config(text="Ganzz. Koeff. werden angestrebt")
            elif method_name == "from_extrema_locations_and_values":
                if aim_int_y_active and "y1_val" in self.param_entries: self.param_row_widgets[1]['hint'].config(text="Ganzzahliger y-Wert") 
                if aim_int_y_active and "y2_val" in self.param_entries: self.param_row_widgets[3]['hint'].config(text="Ganzzahliger y-Wert") 
                if force_int_active and "q1_val" in self.param_entries: 
                    self.param_row_widgets[0]['hint'].config(text="Ganzz. Koeff. werden angestrebt") 
            elif method_name == "from_inflection_point_and_delta_q":
                if aim_int_y_active and "C0" in self.param_entries: self.param_row_widgets[3]['hint'].config(text="Ganzzahlig für ganzz. y-Werte (beeinflusst d)")
                if "K_user" in self.param_entries: self.param_row_widgets[2]['hint'].config(text="K für P'. Vielfaches von 6 für ganzz. Koeff.")
            elif method_name == "from_integer_inflection_real_extrema":
                if "C1_coeff" in self.param_entries: self.param_row_widgets[3]['hint'].config(text="Muss anderes Vorzeichen als 'a' haben für Extrema") 
                if aim_int_y_active and "D0" in self.param_entries: self.param_row_widgets[1]['hint'].config(text="Ganzzahlig für ganzz. y am WP") 
        
        elif poly_type == "Quadratisch":
            method_name = self.quadratic_construction_options.get(construction_key)
            if method_name == "from_quadratic_vertex_and_scale":
                if aim_int_y_active and "k_vertex" in self.param_entries: self.param_row_widgets[1]['hint'].config(text="Ganzzahlig für ganzz. y am Scheitel") 
                if force_int_active and "scale_sq" in self.param_entries: self.param_row_widgets[2]['hint'].config(text="Ganzz. Koeff. wenn h,k,scale_sq ganzz. oder passende Brüche") 
        
        elif poly_type == "Linear":
            method_name = self.linear_construction_options.get(construction_key)
            if method_name == "from_linear_root_and_y_intercept":
                if "root_x_lin" in self.param_entries: self.param_row_widgets[0]['hint'].config(text="Wenn 0, muss y-Abschnitt 0 sein (Steigung dann unbestimmt oder 0)") 
        
        elif poly_type == "Biquadratisch": 
            method_name = self.biquadratic_construction_options.get(construction_key)
            if method_name == "from_outer_extrema_and_central_extremum_biquad":
                if "x_outer_ext" in self.param_entries: self.param_row_widgets[0]['hint'].config(text="x_e > 0")
                if aim_int_y_active and "y_outer_ext" in self.param_entries: self.param_row_widgets[1]['hint'].config(text="Ganzzahliger y-Wert")
                if aim_int_y_active and "y_central_ext" in self.param_entries: self.param_row_widgets[2]['hint'].config(text="Ganzzahliger y-Wert")
                if force_int_active : 
                    hint_text = self.param_row_widgets[0]['hint'].cget("text")
                    self.param_row_widgets[0]['hint'].config(text=(hint_text + " (Ganzz. Koeff. angestrebt)" if hint_text else "(Ganzz. Koeff. angestrebt)"))

            elif method_name == "from_inflection_points_and_central_extremum_biquad":
                if "x_wp" in self.param_entries: self.param_row_widgets[0]['hint'].config(text="x_wp > 0")
                if aim_int_y_active and "y_wp" in self.param_entries: self.param_row_widgets[1]['hint'].config(text="Ganzzahliger y-Wert")
                if aim_int_y_active and "y_central_ext_wp" in self.param_entries: self.param_row_widgets[2]['hint'].config(text="Ganzzahliger y-Wert")
                if force_int_active : 
                    hint_text = self.param_row_widgets[0]['hint'].cget("text")
                    self.param_row_widgets[0]['hint'].config(text=(hint_text + " (Ganzz. Koeff. angestrebt)" if hint_text else "(Ganzz. Koeff. angestrebt)"))
            elif method_name == "direct_coeffs_biquad":
                 if force_int_active and "A_biquad" in self.param_entries : self.param_row_widgets[0]['hint'].config(text="Ganzz. Koeff. angestrebt")
        
        # Hint für x-shift (falls vorhanden)
        if hasattr(self, 'x_shift_specific_param_name') and self.x_shift_specific_param_name in self.param_entries:
            # Finde die Zeile des x_shift-Parameters
            for i in range(self.MAX_PARAM_ROWS):
                if self.param_row_widgets[i]['label'].cget("text") == "X-Verschiebung (s):":
                    self.param_row_widgets[i]['hint'].config(text="Verschiebt Graphen um s")
                    break


    def _update_cubic_construction_params_create_fields(self):
        self.clear_param_input_fields_display() 
        self.input_params_frame.config(text="Parameter für Kubische Konstruktion") 
        construction_key = self.construction_method_var.get()
        method_name = self.cubic_construction_options.get(construction_key)
        current_row = 0
        
        if method_name in ["from_extrema_locations", "from_extrema_locations_and_one_y_target", "from_inflection_point_and_delta_q"]:
            self.flip_extrema_button.pack(side=tk.LEFT, padx=5)
        else:
            self.flip_extrema_button.pack_forget()

        if method_name == "direct_coeffs_cubic": 
            current_row = self._configure_param_entry_row(current_row, "a_cubic_coeff", "Koeffizient a (für x^3):", "1")
            current_row = self._configure_param_entry_row(current_row, "b_cubic_coeff", "Koeffizient b (für x^2):", "0")
            current_row = self._configure_param_entry_row(current_row, "c_cubic_coeff", "Koeffizient c (für x):", "-3")
            current_row = self._configure_param_entry_row(current_row, "d_cubic_coeff", "Koeffizient d (Konstante):", "1")
        elif method_name == "from_saddle_point":
            current_row = self._configure_param_entry_row(current_row, "x0", "x-Sattelpunkt (x0):", "1")
            current_row = self._configure_param_entry_row(current_row, "D0", "y-Sattelpunkt (D0):", "2")
            current_row = self._configure_param_entry_row(current_row, "P_triple_prime_0", "P'''(x0) (6a):", "6")
        elif method_name == "from_extrema_locations":
            current_row = self._configure_param_entry_row(current_row, "q1", "x-Extremum 1 (q1):", "0")
            current_row = self._configure_param_entry_row(current_row, "q2", "x-Extremum 2 (q2):", "2")
            current_row = self._configure_param_entry_row(current_row, "K_user", "Skalierungsfaktor K (für P'):", "1") 
            current_row = self._configure_param_entry_row(current_row, "C0", "Konstante d (y-Versch.):", "0")
        elif method_name == "from_extrema_locations_and_one_y_target": 
            current_row = self._configure_param_entry_row(current_row, "q1_onetarget", "x-Extremum 1 (q1):", "0")
            current_row = self._configure_param_entry_row(current_row, "q2_onetarget", "x-Extremum 2 (q2):", "2")
            current_row = self._configure_param_entry_row(current_row, "y_at_q1_onetarget", "Ziel y-Wert bei q1:", "3")
            current_row = self._configure_param_entry_row(current_row, "K_user_onetarget", "Skalierungsfaktor K (für P'):", "1")
        elif method_name == "from_extrema_locations_and_values": 
            current_row = self._configure_param_entry_row(current_row, "q1_val", "x-Extremum 1 (q1):", "0")
            current_row = self._configure_param_entry_row(current_row, "y1_val", "y-Wert bei q1 (y1):", "0")
            current_row = self._configure_param_entry_row(current_row, "q2_val", "x-Extremum 2 (q2):", "2")
            current_row = self._configure_param_entry_row(current_row, "y2_val", "y-Wert bei q2 (y2):", "4")
        elif method_name == "from_inflection_point_and_delta_q":
            current_row = self._configure_param_entry_row(current_row, "xw", "x-Wendepunkt (xw):", "1")
            current_row = self._configure_param_entry_row(current_row, "delta_q", "Abstand xw zu Extrema (Δq):", "1")
            current_row = self._configure_param_entry_row(current_row, "K_user", "Skalierungsfaktor K (für P'):", "1")
            current_row = self._configure_param_entry_row(current_row, "C0", "Konstante d (y-Versch.):", "0")
        elif method_name == "from_integer_inflection_real_extrema":
            current_row = self._configure_param_entry_row(current_row, "xw", "x-Wendepunkt (xw):", "0")
            current_row = self._configure_param_entry_row(current_row, "D0", "y-Wendepunkt (D0):", "0")
            current_row = self._configure_param_entry_row(current_row, "a_coeff", "Koeff. a (von x^3):", "1")
            current_row = self._configure_param_entry_row(current_row, "C1_coeff", "Koeff. P'(xw):", "-3")
        elif method_name == "from_zero_and_two_integer_roots":
            current_row = self._configure_param_entry_row(current_row, "r2", "Nullstelle r2 (neben x=0):", "1")
            current_row = self._configure_param_entry_row(current_row, "r3", "Nullstelle r3 (neben x=0):", "2")
            current_row = self._configure_param_entry_row(current_row, "scale_factor_S", "Skalierungsfaktor S:", "1")
        
        current_row = self._add_x_shift_parameter_if_applicable(current_row)
        self._apply_param_hints()

    def _update_quadratic_construction_params_create_fields(self):
        self.clear_param_input_fields_display()
        self.input_params_frame.config(text="Parameter für Quadratische Konstruktion")
        construction_key = self.construction_method_var.get()
        method_name = self.quadratic_construction_options.get(construction_key)
        current_row = 0
        self.flip_extrema_button.pack_forget() 

        if method_name == "direct_coeffs_quad": 
            current_row = self._configure_param_entry_row(current_row, "b_quad", "Koeffizient b (für x^2):", "1")
            current_row = self._configure_param_entry_row(current_row, "c_quad", "Koeffizient c (für x):", "0")
            current_row = self._configure_param_entry_row(current_row, "d_quad", "Koeffizient d (Konstante):", "0")
        elif method_name == "from_quadratic_vertex_and_scale":
            current_row = self._configure_param_entry_row(current_row, "h_vertex", "Scheitelpunkt x (h):", "1")
            current_row = self._configure_param_entry_row(current_row, "k_vertex", "Scheitelpunkt y (k):", "2")
            current_row = self._configure_param_entry_row(current_row, "scale_sq", "Skalierungsfaktor (b von x^2):", "1")
        
        self._apply_param_hints()

    def _update_linear_construction_params_create_fields(self):
        self.clear_param_input_fields_display()
        self.input_params_frame.config(text="Parameter für Lineare Konstruktion")
        construction_key = self.construction_method_var.get()
        method_name = self.linear_construction_options.get(construction_key)
        current_row = 0
        self.flip_extrema_button.pack_forget() 

        if method_name == "direct_coeffs_lin": 
            current_row = self._configure_param_entry_row(current_row, "c_lin", "Koeffizient c (Steigung m):", "1")
            current_row = self._configure_param_entry_row(current_row, "d_lin", "Koeffizient d (y-Abschnitt n):", "0")
        elif method_name == "from_linear_root_and_y_intercept":
            current_row = self._configure_param_entry_row(current_row, "root_x_lin", "Nullstelle x:", "1")
            current_row = self._configure_param_entry_row(current_row, "y_intercept_lin", "y-Achsenabschnitt (d):", "2")
        
        self._apply_param_hints()

    def _update_biquadratic_construction_params_create_fields(self): 
        self.clear_param_input_fields_display()
        self.input_params_frame.config(text="Parameter für Biquadratische Konstruktion (Ax^4+Bx^2+D)")
        construction_key = self.construction_method_var.get()
        method_name = self.biquadratic_construction_options.get(construction_key)
        current_row = 0
        self.flip_extrema_button.pack_forget() 

        if method_name == "from_outer_extrema_and_central_extremum_biquad":
            current_row = self._configure_param_entry_row(current_row, "x_outer_ext", "x äußeres Extremum (x_e > 0):", "2")
            current_row = self._configure_param_entry_row(current_row, "y_outer_ext", "y äußeres Extremum (y_e):", "-3")
            current_row = self._configure_param_entry_row(current_row, "y_central_ext", "y zentrales Extremum (y_0 bei x=0):", "1")
        elif method_name == "from_inflection_points_and_central_extremum_biquad":
            current_row = self._configure_param_entry_row(current_row, "x_wp", "x Wendepunkt (x_wp > 0):", "1")
            current_row = self._configure_param_entry_row(current_row, "y_wp", "y Wendepunkt (y_wp):", "0")
            current_row = self._configure_param_entry_row(current_row, "y_central_ext_wp", "y zentrales Extremum (y_0 bei x=0):", "2") 
        elif method_name == "direct_coeffs_biquad":
            current_row = self._configure_param_entry_row(current_row, "A_biquad", "Koeffizient A (für x^4):", "1")
            current_row = self._configure_param_entry_row(current_row, "B_biquad", "Koeffizient B (für x^2):", "-5")
            current_row = self._configure_param_entry_row(current_row, "D_biquad", "Koeffizient D (Konstante):", "4")
        
        current_row = self._add_x_shift_parameter_if_applicable(current_row)
        self._apply_param_hints()

    def get_param_value(self, param_name, param_type_hint=Fraction, is_optional=False, default_if_missing=0):
        try:
            if param_name not in self.param_entries:
                if is_optional: return Fraction(default_if_missing) 
                raise ValueError(f"Parameter '{param_name}' ist für die aktuelle Auswahl nicht verfügbar oder wurde nicht korrekt konfiguriert.")
            val_str = self.param_entries[param_name].get().strip()
            if not val_str: 
                if is_optional: return Fraction(default_if_missing)
                raise ValueError(f"Parameter '{param_name}' darf nicht leer sein.")
            
            f_val = Fraction(val_str)
            if param_type_hint == int:
                if f_val.denominator != 1: raise ValueError(f"Parameter '{param_name}' muss eine ganze Zahl sein, nicht '{val_str}'.")
                return f_val.numerator 
            return f_val 
        except KeyError: 
            if is_optional: return Fraction(default_if_missing)
            raise ValueError(f"Parameter '{param_name}' nicht im UI gefunden (Programmierfehler).")
        except (ValueError, TypeError) as e: 
            if is_optional and not val_str: return Fraction(default_if_missing) 
            original_msg = str(e)
            if param_name not in original_msg:
                 original_msg = f"für Parameter '{param_name}': {original_msg}"
            raise ValueError(f"Ungültiger Wert {original_msg}. Erwartet Zahl oder Bruch (z.B. 1/3 oder 0.5).")

    def generate_and_display(self, use_specific_x_shift=True, called_from_lucky_button=False, parsed_coeffs=None): 
        for key in self.output_fields: 
            self.output_fields[key].config(state=tk.NORMAL); self.output_fields[key].delete(1.0, tk.END); self.output_fields[key].config(state=tk.DISABLED)
        
        base_poly_object = None 
        self.analyzed_poly = None 
        actual_x_shift_used = Fraction(0)

        error_msg_template_cubic = (
            "Erzwungene ganzzahlige Koeffizienten führen dazu, dass der x^3-Term (a) Null wird. "
            "Dies würde den Grad des Polynoms reduzieren. "
            "Bitte Parameter anpassen oder 'Ganzzahlige Koeffizienten' deaktivieren, um ein kubisches Polynom zu erhalten."
        )
        error_msg_template_biquad_A = (
            "Erzwungene ganzzahlige Koeffizienten führen dazu, dass der x^4-Term (A) Null wird. "
            "Dies würde den Grad des Polynoms reduzieren. "
            "Bitte Parameter anpassen oder 'Ganzzahlige Koeffizienten' deaktivieren, um ein biquadratisches Polynom zu erhalten."
        )
        error_msg_template_quad_B = (
            "Erzwungene ganzzahlige Koeffizienten führen dazu, dass der x^2-Term (b) Null wird. "
            "Dies würde den Grad des quadratischen Polynoms reduzieren. "
            "Bitte Parameter anpassen oder 'Ganzzahlige Koeffizienten' deaktivieren."
        )

        try:
            poly_type = self.poly_type_var.get()
            force_int_coeffs_global = self.force_int_coeffs_var.get()
            
            if parsed_coeffs: # Wenn von String-Parser kommend
                actual_x_shift_used = Fraction(0) # Keine zusätzliche Verschiebung
                base_poly_object = Polynomial(
                    parsed_coeffs.get('e',0), 
                    parsed_coeffs.get('a',0), 
                    parsed_coeffs.get('b',0), 
                    parsed_coeffs.get('c',0), 
                    parsed_coeffs.get('d',0)
                )
            else: # Reguläre Generierung über Parameter
                if use_specific_x_shift and poly_type in ["Kubisch", "Biquadratisch"]:
                    actual_x_shift_used = self.get_param_value(self.x_shift_specific_param_name, is_optional=True, default_if_missing=0)
                else: 
                    actual_x_shift_used = Fraction(0)

                # 1. Basis-Polynom erstellen (immer P(x), nicht P(x-s) an dieser Stelle)
                if poly_type == "Kubisch":
                    construction_key, method_name = self.construction_method_var.get(), self.cubic_construction_options.get(self.construction_method_var.get())
                    if not construction_key or not method_name: 
                        if called_from_lucky_button: raise ValueError("Keine Konstruktionsmethode für Kubisch gewählt.")
                        messagebox.showerror("Auswahlfehler", "Bitte Konstruktionsmethode für Kubisch wählen."); return
                    
                    if method_name == "direct_coeffs_cubic":
                        a_val = self.get_param_value("a_cubic_coeff") 
                        b_val = self.get_param_value("b_cubic_coeff")
                        c_val = self.get_param_value("c_cubic_coeff")
                        d_val = self.get_param_value("d_cubic_coeff")
                        if force_int_coeffs_global:
                            if a_val != 0 and a_val.limit_denominator(1) == 0: raise ValueError(error_msg_template_cubic)
                            base_poly_object = Polynomial(0, a_val.limit_denominator(1), b_val.limit_denominator(1), c_val.limit_denominator(1), d_val.limit_denominator(1))
                        else:
                            if a_val == 0 and not called_from_lucky_button: 
                                messagebox.showwarning("Hinweis", "Koeffizient 'a' für x^3 ist Null. Dies ist ein Polynom niedrigeren Grades.")
                            base_poly_object = Polynomial(0,a_val,b_val,c_val,d_val)
                    elif method_name == "from_saddle_point":
                        x0,D0,P_ppp_val = self.get_param_value("x0"), self.get_param_value("D0"), self.get_param_value("P_triple_prime_0") 
                        base_poly_object = from_saddle_point_factory(x0,D0,P_ppp_val,force_int_coeffs_global)
                    elif method_name == "from_extrema_locations":
                        q1,q2,K,C0 = self.get_param_value("q1"),self.get_param_value("q2"),self.get_param_value("K_user"),self.get_param_value("C0")
                        temp_poly = from_extrema_locations_factory(q1,q2,K,C0)
                        if force_int_coeffs_global:
                            if temp_poly.a != 0 and temp_poly.a.limit_denominator(1) == 0: raise ValueError(error_msg_template_cubic)
                            base_poly_object = Polynomial(0, temp_poly.a.limit_denominator(1), temp_poly.b.limit_denominator(1), temp_poly.c.limit_denominator(1), temp_poly.d.limit_denominator(1))
                        else: base_poly_object = temp_poly
                    elif method_name == "from_extrema_locations_and_one_y_target": 
                        q1,q2,y_at_q1,K_user = self.get_param_value("q1_onetarget"),self.get_param_value("q2_onetarget"),self.get_param_value("y_at_q1_onetarget"),self.get_param_value("K_user_onetarget")
                        base_poly_object = from_extrema_locations_and_one_y_target_factory(q1,q2,y_at_q1,K_user,force_int_coeffs_global)
                    elif method_name == "from_extrema_locations_and_values": 
                        q1,y1,q2,y2 = self.get_param_value("q1_val"),self.get_param_value("y1_val"),self.get_param_value("q2_val"),self.get_param_value("y2_val")
                        base_poly_object = from_extrema_locations_and_values_factory(q1,q2,y1,y2,force_int_coeffs_global)
                    elif method_name == "from_inflection_point_and_delta_q":
                        xw,dq,K,C0 = self.get_param_value("xw"),self.get_param_value("delta_q"),self.get_param_value("K_user"),self.get_param_value("C0")
                        temp_poly = from_inflection_point_and_delta_q_factory(xw,dq,K,C0)
                        if force_int_coeffs_global:
                            if temp_poly.a != 0 and temp_poly.a.limit_denominator(1) == 0: raise ValueError(error_msg_template_cubic)
                            base_poly_object = Polynomial(0, temp_poly.a.limit_denominator(1), temp_poly.b.limit_denominator(1), temp_poly.c.limit_denominator(1), temp_poly.d.limit_denominator(1))
                        else: base_poly_object = temp_poly
                    elif method_name == "from_integer_inflection_real_extrema": 
                        xw,D0,a_c,C1_c = self.get_param_value("xw"),self.get_param_value("D0"),self.get_param_value("a_coeff"),self.get_param_value("C1_coeff")
                        temp_poly = from_integer_inflection_real_extrema_factory(xw,D0,a_c,C1_c)
                        if force_int_coeffs_global: 
                             if temp_poly.a != 0 and temp_poly.a.limit_denominator(1) == 0: raise ValueError(error_msg_template_cubic)
                             base_poly_object = Polynomial(0, temp_poly.a.limit_denominator(1), temp_poly.b.limit_denominator(1), temp_poly.c.limit_denominator(1), temp_poly.d.limit_denominator(1))
                        else: base_poly_object = temp_poly 
                    elif method_name == "from_zero_and_two_integer_roots":
                        r2,r3,S = self.get_param_value("r2"),self.get_param_value("r3"),self.get_param_value("scale_factor_S")
                        temp_poly = from_zero_and_two_integer_roots_factory(r2,r3,S)
                        if force_int_coeffs_global:
                             if temp_poly.a != 0 and temp_poly.a.limit_denominator(1) == 0: raise ValueError(error_msg_template_cubic)
                             base_poly_object = Polynomial(0, temp_poly.a.limit_denominator(1), temp_poly.b.limit_denominator(1), temp_poly.c.limit_denominator(1), temp_poly.d.limit_denominator(1))
                        else: base_poly_object = temp_poly
                
                elif poly_type == "Quadratisch":
                    construction_key_quad, method_name_quad = self.construction_method_var.get(), self.quadratic_construction_options.get(self.construction_method_var.get())
                    if not construction_key_quad or not method_name_quad: 
                        if called_from_lucky_button: raise ValueError("Keine Konstruktionsmethode für Quadratisch gewählt.")
                        messagebox.showerror("Auswahlfehler", "Bitte Konstruktionsmethode für Quadratisch wählen."); return
                    
                    if method_name_quad == "from_quadratic_vertex_and_scale":
                        h,k,scale_sq = self.get_param_value("h_vertex"),self.get_param_value("k_vertex"),self.get_param_value("scale_sq")
                        base_poly_object = from_quadratic_vertex_and_scale_factory(h,k,scale_sq,force_int_coeffs_global) 
                    elif method_name_quad == "direct_coeffs_quad":
                        b_param,c_param,d_param = self.get_param_value("b_quad"),self.get_param_value("c_quad"),self.get_param_value("d_quad")
                        if force_int_coeffs_global:
                            original_b_quad = b_param
                            b_int_quad = original_b_quad.limit_denominator(1)
                            if original_b_quad != 0 and b_int_quad == 0:
                                 raise ValueError(error_msg_template_quad_B)
                            base_poly_object = Polynomial(0,0,b_int_quad,c_param.limit_denominator(1),d_param.limit_denominator(1))
                        else:
                            if b_param == 0 and not called_from_lucky_button: messagebox.showwarning("Hinweis", "Koeffizient 'b' für x^2 ist Null. Dies ist ein Polynom niedrigeren Grades.")
                            base_poly_object = Polynomial(0,0,b_param,c_param,d_param) 
                
                elif poly_type == "Linear":
                    construction_key_lin, method_name_lin = self.construction_method_var.get(), self.linear_construction_options.get(self.construction_method_var.get())
                    if not construction_key_lin or not method_name_lin: 
                        if called_from_lucky_button: raise ValueError("Keine Konstruktionsmethode für Linear gewählt.")
                        messagebox.showerror("Auswahlfehler", "Bitte Konstruktionsmethode für Linear wählen."); return
                    if method_name_lin == "from_linear_root_and_y_intercept":
                        root_x,y_int = self.get_param_value("root_x_lin"),self.get_param_value("y_intercept_lin")
                        base_poly_object = from_linear_root_and_y_intercept_factory(root_x,y_int) 
                    elif method_name_lin == "direct_coeffs_lin":
                        c,d = self.get_param_value("c_lin"),self.get_param_value("d_lin")
                        if force_int_coeffs_global: 
                             base_poly_object = Polynomial(0,0,0,c.limit_denominator(1),d.limit_denominator(1))
                        else:
                             if c == 0 and not called_from_lucky_button: messagebox.showwarning("Hinweis", "Koeffizient 'c' (Steigung) ist Null. Dies ist eine konstante Funktion.")
                             base_poly_object = Polynomial(0,0,0,c,d) 

                elif poly_type == "Biquadratisch": 
                    construction_key_biq, method_name_biq = self.construction_method_var.get(), self.biquadratic_construction_options.get(self.construction_method_var.get())
                    if not construction_key_biq or not method_name_biq: 
                        if called_from_lucky_button: raise ValueError("Keine Konstruktionsmethode für Biquadratisch gewählt.")
                        messagebox.showerror("Auswahlfehler", "Bitte Konstruktionsmethode für Biquadratisch wählen."); return

                    if method_name_biq == "from_outer_extrema_and_central_extremum_biquad":
                        xe, ye, y0 = self.get_param_value("x_outer_ext"), self.get_param_value("y_outer_ext"), self.get_param_value("y_central_ext")
                        base_poly_object = from_outer_extrema_and_central_extremum_biquad_factory(xe, ye, y0, force_int_coeffs_global) 
                    elif method_name_biq == "from_inflection_points_and_central_extremum_biquad":
                        xwp, ywp, y0_wp = self.get_param_value("x_wp"), self.get_param_value("y_wp"), self.get_param_value("y_central_ext_wp")
                        base_poly_object = from_inflection_points_and_central_extremum_biquad_factory(xwp, ywp, y0_wp, force_int_coeffs_global) 
                    elif method_name_biq == "direct_coeffs_biquad":
                        A_param,B_param,D_param = self.get_param_value("A_biquad"), self.get_param_value("B_biquad"), self.get_param_value("D_biquad")
                        if force_int_coeffs_global:
                             original_A_biq = A_param
                             A_int_biq = original_A_biq.limit_denominator(1)
                             if original_A_biq != 0 and A_int_biq == 0:
                                 raise ValueError(error_msg_template_biquad_A)
                             base_poly_object = direct_coeffs_biquad_factory(A_int_biq, B_param.limit_denominator(1), D_param.limit_denominator(1))
                        else:
                             if A_param == 0 and not called_from_lucky_button: messagebox.showwarning("Hinweis", "Koeffizient 'A' für x^4 ist Null. Dies ist ein Polynom niedrigeren Grades.")
                             base_poly_object = direct_coeffs_biquad_factory(A_param,B_param,D_param)

            # 2. Wende die Verschiebung an, um das finale Polynom zu erhalten (nur wenn nicht von String-Parser)
            if base_poly_object:
                if actual_x_shift_used != 0 and not parsed_coeffs: # Nur verschieben, wenn nicht schon durch Parser definiert
                    self.analyzed_poly = base_poly_object.expand_and_shift(actual_x_shift_used)
                else:
                    self.analyzed_poly = base_poly_object # Ist bereits das finale Polynom
            else:
                 if not called_from_lucky_button: messagebox.showerror("Fehler", "Basis-Polynom konnte nicht erstellt werden.")
                 else: raise ValueError("Basis-Polynom konnte nicht erstellt werden (Lucky Button).")
                 return

            if self.analyzed_poly: 
                self.display_poly_info(self.analyzed_poly) 
                self.export_discussion_button.config(state=tk.NORMAL) 
                self.update_status("Polynom generiert und analysiert.")
            else: 
                if not called_from_lucky_button: messagebox.showerror("Fehler", "Polynom konnte nicht generiert werden (nach Verschiebung).")
                else: raise ValueError("Polynom konnte nicht generiert werden nach Verschiebung (Lucky Button).")
                self.export_discussion_button.config(state=tk.DISABLED)
        except ValueError as e: 
            if called_from_lucky_button: raise e 
            messagebox.showerror("Eingabefehler oder ungültige Konstruktion", str(e))
            self.export_discussion_button.config(state=tk.DISABLED)
            self.update_status(f"Fehler: {str(e)[:50]}...") # Kurze Fehlermeldung in Statusleiste
        except Exception as e: 
            import traceback
            messagebox.showerror("Unerwarteter Fehler", f"Fehler: {str(e)}\n\nTraceback:\n{traceback.format_exc()}")
            self.export_discussion_button.config(state=tk.DISABLED)
            self.update_status("Unerwarteter Fehler.")

    def _format_number_list_for_display(self, data_list, precision=3):
        if data_list is None or (isinstance(data_list, list) and not data_list):
            return "Keine"
        if isinstance(data_list, list):
            if data_list and isinstance(data_list[0], str): 
                return data_list[0]
            try:
                return ", ".join([f"{float(x):.{precision}f}" for x in data_list if isinstance(x, (int, float, Fraction))])
            except (ValueError, TypeError):
                return "Formatierungsfehler" 
        try: 
            return f"{float(data_list):.{precision}f}"
        except (ValueError, TypeError):
            return "Formatierungsfehler"


    def display_poly_info(self, poly_to_display): 
        
        p_str = str(poly_to_display) 
        
        deriv1_coeffs = poly_to_display.derivative1_coeffs() 
        deriv2_coeffs = poly_to_display.derivative2_coeffs() 

        p_prime_obj = Polynomial.from_coeffs(deriv1_coeffs[1:]) 
        p_double_prime_obj = Polynomial.from_coeffs(deriv2_coeffs[2:])

        p_prime_str = str(p_prime_obj)
        p_double_prime_str = str(p_double_prime_obj)

        roots_found = poly_to_display.find_roots() 
        crit_x_float_list = poly_to_display.find_critical_points_x() 
        infl_x_float_list = poly_to_display.find_inflection_points_x() 

        roots_str = self._format_number_list_for_display(roots_found)
        crit_x_str = self._format_number_list_for_display(crit_x_float_list)
        
        crit_y_str = "Keine"
        if isinstance(crit_x_float_list, list) and crit_x_float_list and not isinstance(crit_x_float_list[0], str):
            crit_y_values = [poly_to_display.evaluate(x) for x in crit_x_float_list] 
            crit_y_str = self._format_number_list_for_display(crit_y_values)
        elif crit_x_str != "Keine" and crit_x_str != "Alle reellen Zahlen" and crit_x_str != "Formatierungsfehler": 
            try:
                crit_y_str = self._format_number_list_for_display([poly_to_display.evaluate(float(crit_x_str))])
            except: pass 

        infl_x_str = self._format_number_list_for_display(infl_x_float_list)
        infl_y_str = "Keine"
        if isinstance(infl_x_float_list, list) and infl_x_float_list and not isinstance(infl_x_float_list[0], str):
            infl_y_values = [poly_to_display.evaluate(x) for x in infl_x_float_list] 
            infl_y_str = self._format_number_list_for_display(infl_y_values)
        elif infl_x_str != "Keine" and infl_x_str != "Alle reellen Zahlen" and infl_x_str != "Formatierungsfehler":
            try:
                infl_y_str = self._format_number_list_for_display([poly_to_display.evaluate(float(infl_x_str))])
            except: pass


        saddle_info_str = "Kein Sattelpunkt"
        if isinstance(crit_x_float_list, list) and isinstance(infl_x_float_list, list) and \
           not (crit_x_float_list and isinstance(crit_x_float_list[0], str)) and \
           not (infl_x_float_list and isinstance(infl_x_float_list[0], str)):
            saddle_points_x = []
            for infl_x in infl_x_float_list:
                if any(abs(infl_x - crit_x) < 1e-7 for crit_x in crit_x_float_list):
                    saddle_points_x.append(infl_x) 
            if saddle_points_x:
                saddle_info_str = "Ja, bei x = " + self._format_number_list_for_display(saddle_points_x)
        elif not isinstance(infl_x_float_list, list) or not infl_x_float_list or (isinstance(infl_x_float_list[0], str) and infl_x_float_list[0] == "Keine"):
            saddle_info_str = "Kein Sattelpunkt (da keine WP)"


        outputs = {"P(x)": p_str, "P'(x)": p_prime_str, "P''(x)": p_double_prime_str, 
                   "roots": roots_str, "crit_x": crit_x_str, "crit_y": crit_y_str, 
                   "infl_x": infl_x_str, "infl_y": infl_y_str, "saddle_info": saddle_info_str}
        
        for key, value in outputs.items():
            self.output_fields[key].config(state=tk.NORMAL); self.output_fields[key].delete(1.0, tk.END)
            self.output_fields[key].insert(tk.END, value); self.output_fields[key].config(state=tk.DISABLED)

        # Plot aktualisieren
        self.plot_ax.clear() 
        
        initial_view_center_x = 0 
        all_important_x_for_view = [] 
        if roots_found and isinstance(roots_found, list) and not (roots_found and isinstance(roots_found[0], str)):
            all_important_x_for_view.extend([float(r) for r in roots_found if isinstance(r, (int, float, Fraction))])
        if crit_x_float_list and isinstance(crit_x_float_list, list) and not (crit_x_float_list and isinstance(crit_x_float_list[0], str)):
             all_important_x_for_view.extend(crit_x_float_list)
        if infl_x_float_list and isinstance(infl_x_float_list, list) and not (infl_x_float_list and isinstance(infl_x_float_list[0], str)):
            all_important_x_for_view.extend(infl_x_float_list)
        
        if all_important_x_for_view:
            initial_view_center_x = (min(all_important_x_for_view) + max(all_important_x_for_view)) / 2
            min_x_val_s, max_x_val_s = min(all_important_x_for_view), max(all_important_x_for_view)
            delta_x_s = max_x_val_s - min_x_val_s
            padding_s = max(delta_x_s * 0.5, 3.0) 
            if self.poly_type_var.get() == "Kubisch" and self.cubic_construction_options.get(self.construction_method_var.get()) == "from_saddle_point":
                padding_s = max(padding_s, 4.0) 
            x_lim_min, x_lim_max = min_x_val_s - padding_s, max_x_val_s + padding_s
        else: 
            x_lim_min, x_lim_max = initial_view_center_x - 5, initial_view_center_x + 5 

        calc_x_min = initial_view_center_x - 50 
        calc_x_max = initial_view_center_x + 50
        x_vals_for_plot_calc = np.linspace(calc_x_min, calc_x_max, 800) 
        
        y_vals_calc = poly_to_display.evaluate(x_vals_for_plot_calc) 
        y_vals_calc_float = np.array(y_vals_calc, dtype=float)

        plot_label_str = "P(x)" 
        self.plot_ax.plot(x_vals_for_plot_calc, y_vals_calc_float, label=plot_label_str, color='blue')
        
        if roots_found and isinstance(roots_found, list) and not (roots_found and isinstance(roots_found[0], str)):
            roots_plot = [float(r) for r in roots_found if isinstance(r, (int, float, Fraction))]
            if roots_plot:
                 self.plot_ax.scatter(roots_plot, [0]*len(roots_plot), color='black', s=50, zorder=5, label="Nullstellen")

        if crit_x_float_list and isinstance(crit_x_float_list, list) and not (crit_x_float_list and isinstance(crit_x_float_list[0], str)):
            saddle_points_x_for_plot = []
            if 'saddle_points_x' in locals() and saddle_points_x: 
                 saddle_points_x_for_plot = saddle_points_x

            non_saddle_crit_x = [x for x in crit_x_float_list if x not in saddle_points_x_for_plot]
            if non_saddle_crit_x:
                non_saddle_crit_y_values = [float(poly_to_display.evaluate(c_val)) for c_val in non_saddle_crit_x]
                if len(non_saddle_crit_x) == len(non_saddle_crit_y_values):
                    self.plot_ax.scatter(non_saddle_crit_x, non_saddle_crit_y_values, color='red', zorder=5, s=70, label="Extrema", edgecolors='black')
        
        if infl_x_float_list and isinstance(infl_x_float_list, list) and not (infl_x_float_list and isinstance(infl_x_float_list[0], str)):
            for idx, infl_x_val in enumerate(infl_x_float_list):
                is_saddle_plot = 'saddle_points_x' in locals() and saddle_points_x and any(abs(infl_x_val - spx_s) < 1e-7 for spx_s in saddle_points_x)
                color_infl = 'purple' if is_saddle_plot else 'green'
                label_infl_plot = "Sattelpunkt" if is_saddle_plot else "Wendepunkt"
                marker_infl = '*' if is_saddle_plot else 'P' 
                size_infl = 200 if is_saddle_plot else 150
                current_labels_plot = [c.get_label() for c in self.plot_ax.collections] + [c.get_label() for c in self.plot_ax.lines]
                y_val_for_infl = float(poly_to_display.evaluate(infl_x_val))
                if label_infl_plot not in current_labels_plot:
                    self.plot_ax.scatter([infl_x_val], [y_val_for_infl], 
                                     color=color_infl, zorder=6, s=size_infl, marker=marker_infl, 
                                     label=label_infl_plot, edgecolors='black', linewidth=1)
                else: 
                     self.plot_ax.scatter([infl_x_val], [y_val_for_infl], 
                                     color=color_infl, zorder=6, s=size_infl, marker=marker_infl, 
                                     edgecolors='black', linewidth=1)

        self.plot_ax.set_title(f"Graph von {p_str}") 
        self.plot_ax.set_xlabel("x")
        self.plot_ax.set_ylabel("P(x)") 
        self.plot_ax.grid(True, linestyle=':', alpha=0.7)
        self.plot_ax.axhline(0, color='black', linewidth=0.5)
        self.plot_ax.axvline(0, color='black', linewidth=0.5) 
        
        self.plot_ax.set_xlim(x_lim_min, x_lim_max)
        if y_vals_calc_float.size > 0 and np.all(np.isfinite(y_vals_calc_float)) : 
            visible_y_mask = (x_vals_for_plot_calc >= x_lim_min) & (x_vals_for_plot_calc <= x_lim_max)
            if np.any(visible_y_mask):
                y_vals_in_view = y_vals_calc_float[visible_y_mask]
                if y_vals_in_view.size > 0:
                    min_y_plot, max_y_plot = np.min(y_vals_in_view), np.max(y_vals_in_view)
                    y_range_plot = max_y_plot - min_y_plot
                    padding_y = y_range_plot * 0.1 if abs(y_range_plot) > 1e-6 else 1.0
                    if padding_y == 0 : padding_y = 1.0 
                    self.plot_ax.set_ylim(min_y_plot - padding_y, max_y_plot + padding_y)
                else: 
                    self.plot_ax.set_ylim(-10, 10)
            else: 
                 self.plot_ax.set_ylim(-10, 10)
        else: 
            self.plot_ax.set_ylim(-10, 10)

        self.plot_ax.legend(loc='best')
        self.canvas.draw()

    def export_curve_discussion(self):
        if not self.analyzed_poly:
            messagebox.showwarning("Kein Polynom", "Bitte zuerst ein Polynom generieren und analysieren.")
            return

        features = PolynomialFeatures(self.analyzed_poly)
        discussion_text = features.get_full_discussion()
        ExportDialog(self.master, "Kurvendiskussion Export", discussion_text)
        self.update_status("Kurvendiskussion exportiert.")
        
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
                    if min_val == 0 and max_val == 0: num = random.choice([-1, 1]) 
                    elif min_val == 0: num_min_actual = 1
                    elif max_val == 0: num_max_actual = -1
                num = random.randint(num_min_actual, num_max_actual)
                if non_zero and num == 0: 
                    possible_nums = [i for i in range(num_min_actual, num_max_actual + 1) if i != 0]
                    if not possible_nums: num = 1 if num_max_actual >=1 else -1 
                    else: num = random.choice(possible_nums)
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
        self.poly_string_input_var.set("") # String-Eingabe leeren
        poly_type = self.poly_type_var.get()
        current_construction_options = {}
        if poly_type == "Kubisch": current_construction_options = self.cubic_construction_options
        elif poly_type == "Biquadratisch": current_construction_options = self.biquadratic_construction_options
        elif poly_type == "Quadratisch": current_construction_options = self.quadratic_construction_options
        elif poly_type == "Linear": current_construction_options = self.linear_construction_options
        
        if not current_construction_options:
            self.generate_and_display(use_specific_x_shift=False, called_from_lucky_button=True) 
            return

        max_retries = 20 
        for attempt in range(max_retries):
            try:
                if hasattr(self, 'x_shift_specific_param_name') and self.x_shift_specific_param_name in self.param_entries:
                    self.param_entries[self.x_shift_specific_param_name].set("0")

                chosen_method_key = random.choice(list(current_construction_options.keys()))
                self.construction_method_var.set(chosen_method_key) 
                self._update_construction_params_display() 

                if poly_type != "Linear":
                    self.force_int_coeffs_var.set(random.choice([True, False]))
                    self.aim_int_y_var.set(random.choice([True, False]))
                self._apply_param_hints() 
                
                x_coord_range = (-3, 3); y_coord_range = (-5, 5); scale_val_range = (-3, 3) 
                delta_range = (1, 2); p_ppp_forced_int_choices = [-18, -12, -6, 6, 12, 18]
                biquad_A_range = (-2,2); biquad_BD_range = (-5,5)

                for param_name, entry_var_tk in self.param_entries.items():
                    if hasattr(self, 'x_shift_specific_param_name') and param_name == self.x_shift_specific_param_name: 
                        entry_var_tk.set("0")
                        continue

                    random_val_str = "1" 
                    is_int_only = False; non_zero_needed = False
                    current_min, current_max = x_coord_range 
                    if poly_type == "Kubisch":
                        if param_name in ["x0", "q1", "q2", "q1_val", "q2_val", "xw", "r2", "r3"]: is_int_only = True 
                        elif param_name in ["D0", "C0", "y_at_q1_onetarget", "y1_val", "y2_val"]:
                            is_int_only = self.aim_int_y_var.get(); current_min, current_max = y_coord_range
                        elif param_name in ["K_user", "a_coeff", "scale_factor_S", "K_user_onetarget","a_cubic_coeff"]: # a_cubic_coeff hinzugefügt
                            non_zero_needed = True; current_min, current_max = scale_val_range
                            if param_name == "K_user" and self.force_int_coeffs_var.get(): 
                                random_val_str = str(random.choice([-18,-12,-6,6,12,18]))
                                entry_var_tk.set(random_val_str); continue
                        elif param_name == "P_triple_prime_0":
                            if self.force_int_coeffs_var.get(): random_val_str = str(random.choice(p_ppp_forced_int_choices))
                            else: non_zero_needed = True; current_min, current_max = scale_val_range
                            entry_var_tk.set(random_val_str); continue
                        elif param_name == "delta_q": is_int_only = True; non_zero_needed = True; current_min, current_max = delta_range
                        elif param_name in ["C1_coeff", "b_cubic_coeff", "c_cubic_coeff", "d_cubic_coeff"]: # b,c,d für kubisch hinzugefügt
                             current_min, current_max = scale_val_range # Können auch 0 sein
                    elif poly_type == "Quadratisch":
                        if param_name in ["h_vertex"]: is_int_only = True
                        elif param_name in ["k_vertex"]: is_int_only = self.aim_int_y_var.get(); current_min, current_max = y_coord_range
                        elif param_name in ["scale_sq", "b_quad"]: non_zero_needed = True; current_min, current_max = scale_val_range
                        elif param_name in ["c_quad", "d_quad"]: current_min, current_max = y_coord_range 
                    elif poly_type == "Linear":
                         if param_name in ["root_x_lin"]: is_int_only = True 
                         elif param_name in ["y_intercept_lin", "d_lin"]: current_min, current_max = y_coord_range
                         elif param_name in ["c_lin"]: non_zero_needed = False; current_min, current_max = scale_val_range 
                    elif poly_type == "Biquadratisch": 
                        if param_name in ["x_outer_ext", "x_wp"]:
                            is_int_only = True; non_zero_needed = True; current_min, current_max = (1,3) 
                        elif param_name in ["y_outer_ext", "y_central_ext", "y_wp", "y_central_ext_wp", "D_biquad"]:
                            is_int_only = self.aim_int_y_var.get(); current_min, current_max = biquad_BD_range
                        elif param_name in ["A_biquad"]:
                            non_zero_needed = True; current_min, current_max = biquad_A_range
                        elif param_name in ["B_biquad"]:
                             current_min, current_max = biquad_BD_range 
                    if not ((poly_type == "Kubisch" and param_name == "P_triple_prime_0") or \
                            (poly_type == "Kubisch" and param_name == "K_user" and self.force_int_coeffs_var.get())):
                        random_val_str = self._generate_random_value(is_int_only, non_zero_needed, current_min, current_max)
                    entry_var_tk.set(random_val_str)

                method_key_lucky = self.construction_method_var.get() 
                if poly_type == "Kubisch":
                    active_method_name = self.cubic_construction_options.get(method_key_lucky)
                    if active_method_name in ["from_extrema_locations", "from_extrema_locations_and_values", "from_extrema_locations_and_one_y_target"]:
                        q1_var_name = "q1" if active_method_name == "from_extrema_locations" else ("q1_val" if active_method_name == "from_extrema_locations_and_values" else "q1_onetarget")
                        q2_var_name = "q2" if active_method_name == "from_extrema_locations" else ("q2_val" if active_method_name == "from_extrema_locations_and_values" else "q2_onetarget")
                        q1_var, q2_var = self.param_entries.get(q1_var_name), self.param_entries.get(q2_var_name)
                        if q1_var and q2_var and q1_var.get() == q2_var.get():
                            q2_var.set(self._generate_random_value(True, False, x_coord_range[0], x_coord_range[1], force_different_from=Fraction(q1_var.get())))
                    elif active_method_name == "from_integer_inflection_real_extrema":
                        a_coeff_var, c1_coeff_var = self.param_entries.get("a_coeff"), self.param_entries.get("C1_coeff")
                        if a_coeff_var and c1_coeff_var:
                            a_val, c1_val = Fraction(a_coeff_var.get()), Fraction(c1_coeff_var.get())
                            if c1_val != 0 and (a_val * c1_val >=0): 
                                c1_coeff_var.set(str(-c1_val if c1_val !=0 else Fraction(random.choice([-1,1])) * abs(a_val))) 
                            elif c1_val == 0 and a_val == 0 : a_coeff_var.set("1") 
                elif poly_type == "Linear":
                    if self.linear_construction_options.get(method_key_lucky) == "from_linear_root_and_y_intercept":
                        root_x_var, y_int_var = self.param_entries.get("root_x_lin"), self.param_entries.get("y_intercept_lin")
                        if root_x_var and y_int_var and Fraction(root_x_var.get()) == 0 and Fraction(y_int_var.get()) != 0:
                            y_int_var.set("0") 
                        elif root_x_var and y_int_var and Fraction(root_x_var.get()) == 0 and Fraction(y_int_var.get()) == 0:
                             root_x_var.set("1") 
                elif poly_type == "Biquadratisch": 
                    active_method_name = self.biquadratic_construction_options.get(method_key_lucky)
                    if active_method_name == "from_outer_extrema_and_central_extremum_biquad":
                        xe_var, ye_var, y0_var = self.param_entries.get("x_outer_ext"), self.param_entries.get("y_outer_ext"), self.param_entries.get("y_central_ext")
                        if xe_var and Fraction(xe_var.get()) == 0: xe_var.set("1") 
                        if ye_var and y0_var and Fraction(ye_var.get()) == Fraction(y0_var.get()): 
                            y0_var.set(str(Fraction(ye_var.get()) + random.choice([-1,1])))
                    elif active_method_name == "from_inflection_points_and_central_extremum_biquad":
                        xwp_var,ywp_var,y0wp_var = self.param_entries.get("x_wp"),self.param_entries.get("y_wp"),self.param_entries.get("y_central_ext_wp")
                        if xwp_var and Fraction(xwp_var.get()) == 0: xwp_var.set("1") 
                        if ywp_var and y0wp_var and Fraction(ywp_var.get()) == Fraction(y0wp_var.get()): 
                            y0wp_var.set(str(Fraction(ywp_var.get()) + random.choice([-1,1])))
                self._apply_param_hints() 
                self.generate_and_display(use_specific_x_shift=False, called_from_lucky_button=True) 
                return 
            except ValueError as e: 
                print(f"I Feel Lucky - Versuch {attempt+1}/{max_retries} fehlgeschlagen (ValueError): {e}")
                if attempt == max_retries - 1:
                    messagebox.showwarning("Zufallsgenerator", 
                                           f"Konnte nach {max_retries} Versuchen kein gültiges Zufallspolynom erzeugen. "
                                           "Möglicherweise sind die internen Zufallsbedingungen zu restriktiv. "
                                           "Bitte Parameter manuell anpassen oder erneut versuchen.")
            except Exception as e_gen: 
                 print(f"I Feel Lucky - Allgemeiner Fehler bei Versuch {attempt+1}: {e_gen}")
                 if attempt == max_retries - 1: 
                    messagebox.showerror("Zufallsgenerator Fehler", f"Unerwarteter Fehler im Zufallsgenerator nach {max_retries} Versuchen: {e_gen}")


    def flip_extrema_type(self):
        poly_type = self.poly_type_var.get()
        if poly_type != "Kubisch":
            messagebox.showinfo("Info", "Extrema-Typ umkehren ist nur für kubische Funktionen mit K_user relevant.")
            return
        construction_key = self.construction_method_var.get()
        method_name = self.cubic_construction_options.get(construction_key)
        param_to_flip = None
        if method_name == "from_extrema_locations": param_to_flip = "K_user"
        elif method_name == "from_extrema_locations_and_one_y_target": param_to_flip = "K_user_onetarget"
        elif method_name == "from_inflection_point_and_delta_q": param_to_flip = "K_user"
        if param_to_flip and param_to_flip in self.param_entries:
            try:
                current_val_str = self.param_entries[param_to_flip].get()
                current_frac = Fraction(current_val_str)
                if current_frac != 0: 
                    self.param_entries[param_to_flip].set(str(-current_frac))
                    self.generate_and_display()
                else:
                    messagebox.showinfo("Info", f"{param_to_flip} ist Null und kann nicht einfach umgekehrt werden, um den Typ zu ändern.")
            except ValueError:
                messagebox.showerror("Fehler", f"Ungültiger Wert im Feld {param_to_flip}.")
        else:
            messagebox.showinfo("Info", "Extrema-Typ umkehren ist für die aktuelle kubische Konstruktionsmethode nicht anwendbar oder der Parameter K fehlt.")

    def parse_term(self, term_str_in):
        term_str = term_str_in.strip()
        
        if not term_str:
            return Fraction(0), 0 # Leerer Term

        # Vorzeichen extrahieren
        sign = 1
        if term_str.startswith('+'):
            term_str = term_str[1:].strip()
        elif term_str.startswith('-'):
            sign = -1
            term_str = term_str[1:].strip()

        if not term_str: # Nur ein Vorzeichen
             raise ValueError("Isolierter Vorzeichen-Term")


        coeff_str = ""
        power_val = 0
        
        if 'x' in term_str:
            parts = term_str.split('x', 1)
            coeff_part = parts[0].strip()
            
            if not coeff_part: # x, -x, +x
                coeff_val = Fraction(sign)
            else:
                try:
                    coeff_val = Fraction(coeff_part) * sign
                except ValueError:
                    raise ValueError(f"Ungültiger Koeffizient: '{coeff_part}'")

            if '^' in parts[1]:
                power_part_str = parts[1].split('^', 1)[1].strip()
                try:
                    power_val = int(power_part_str)
                    if power_val < 0 or power_val > 4 : raise ValueError() # Max Grad 4
                except ValueError:
                    raise ValueError(f"Ungültige Potenz: '^{power_part_str}'")
            elif parts[1].strip() == "": # Term endet mit x
                power_val = 1
            else: # Etwas nach x aber kein ^, z.B. "2xy"
                raise ValueError(f"Unerwartete Zeichen nach 'x': '{parts[1]}'")
        else: # Konstanter Term
            try:
                coeff_val = Fraction(term_str) * sign
                power_val = 0
            except ValueError:
                raise ValueError(f"Ungültiger konstanter Term: '{term_str}'")
                
        return coeff_val, power_val

    def parse_polynomial_string(self, poly_string_raw):
        coeffs = {'e': Fraction(0), 'a': Fraction(0), 'b': Fraction(0), 'c': Fraction(0), 'd': Fraction(0)}
        
        if not poly_string_raw.strip(): # Leere Eingabe
            return coeffs # Alle Koeffizienten bleiben 0

        # Ersetze Leerzeichen um Operatoren herum, aber nicht innerhalb von Brüchen
        poly_string = poly_string_raw.replace(" ", "") 
        poly_string = poly_string.replace("-", "+-") 
        
        if poly_string.startswith("+-"): 
            poly_string = poly_string[1:]
        elif poly_string.startswith("+"): # Falls der User "+2x" etc. eingibt
            poly_string = poly_string[1:]

        if not poly_string: # Wenn nach Normalisierung leer (z.B. nur "-")
            return coeffs

        terms = [term for term in poly_string.split('+') if term] 

        for term_str in terms:
            try:
                coeff, power = self.parse_term(term_str) # parse_term ist jetzt eine Methode der Klasse
                if power == 4: coeffs['e'] += coeff
                elif power == 3: coeffs['a'] += coeff
                elif power == 2: coeffs['b'] += coeff
                elif power == 1: coeffs['c'] += coeff
                elif power == 0: coeffs['d'] += coeff
                else: # Sollte durch parse_term abgefangen werden
                    raise ValueError(f"Potenz {power} wird nicht unterstützt im Term '{term_str}'.")
            except ValueError as e:
                raise ValueError(f"Fehler beim Parsen des Terms '{term_str}': {e}")
        return coeffs

    def apply_parsed_polynomial(self):
        poly_string = self.poly_string_input_var.get()
        if not poly_string.strip():
            messagebox.showinfo("Info", "Bitte geben Sie eine Polynomfunktion ein.")
            return

        try:
            coeffs = self.parse_polynomial_string(poly_string)

            # Polynomtyp bestimmen
            determined_type = None
            target_construction_method = None
            params_to_set = {}

            if coeffs['e'] != 0: # Potenz 4 vorhanden
                if coeffs['a'] == 0 and coeffs['c'] == 0: # Biquadratisch
                    determined_type = "Biquadratisch"
                    target_construction_method = "direct_coeffs_biquad"
                    params_to_set = {"A_biquad": coeffs['e'], "B_biquad": coeffs['b'], "D_biquad": coeffs['d']}
                else:
                    messagebox.showerror("Fehler", "Allgemeine Polynome 4. Grades werden nicht direkt unterstützt. Bitte als biquadratisches Polynom (nur x^4, x^2, Konstante) oder Polynom niedrigeren Grades eingeben.")
                    return
            elif coeffs['a'] != 0: # Höchste Potenz ist 3
                determined_type = "Kubisch"
                target_construction_method = "direct_coeffs_cubic"
                params_to_set = {"a_cubic_coeff": coeffs['a'], "b_cubic_coeff": coeffs['b'], "c_cubic_coeff": coeffs['c'], "d_cubic_coeff": coeffs['d']}
            elif coeffs['b'] != 0: # Höchste Potenz ist 2
                determined_type = "Quadratisch"
                target_construction_method = "direct_coeffs_quad"
                params_to_set = {"b_quad": coeffs['b'], "c_quad": coeffs['c'], "d_quad": coeffs['d']}
            elif coeffs['c'] != 0 or coeffs['d'] != 0 : # Höchste Potenz ist 1 oder 0 (und nicht alles 0)
                determined_type = "Linear"
                target_construction_method = "direct_coeffs_lin"
                params_to_set = {"c_lin": coeffs['c'], "d_lin": coeffs['d']}
            else: # Nullpolynom
                determined_type = "Linear" # Behandle als lineare Funktion P(x)=0
                target_construction_method = "direct_coeffs_lin"
                params_to_set = {"c_lin": 0, "d_lin": 0}
                messagebox.showinfo("Info", "Das eingegebene Polynom ist das Nullpolynom P(x) = 0.")


            if determined_type:
                self.poly_type_var.set(determined_type)
                self.update_ui_for_poly_type() # Aktualisiert die Konstruktionsmethoden-Liste

                # Finde den exakten Schlüssel für die "Aus Koeffizienten"-Methode
                actual_method_key = None
                current_options = {}
                if determined_type == "Kubisch": current_options = self.cubic_construction_options
                elif determined_type == "Biquadratisch": current_options = self.biquadratic_construction_options
                elif determined_type == "Quadratisch": current_options = self.quadratic_construction_options
                elif determined_type == "Linear": current_options = self.linear_construction_options
                
                for key, val_method_name in current_options.items():
                    if val_method_name == target_construction_method:
                        actual_method_key = key
                        break
                
                if actual_method_key:
                    self.construction_method_var.set(actual_method_key)
                    self._update_construction_params_display() # Zeigt die korrekten Parameterfelder an

                    # Parameterfelder füllen
                    for param_key, value_to_set in params_to_set.items():
                        if param_key in self.param_entries:
                            self.param_entries[param_key].set(str(value_to_set))
                        else:
                            print(f"Warnung: Parameter-Schlüssel '{param_key}' nicht in self.param_entries gefunden beim Parsen.")
                    
                    # X-Shift auf 0 setzen für geparste Polynome
                    if hasattr(self, 'x_shift_specific_param_name') and self.x_shift_specific_param_name in self.param_entries:
                         self.param_entries[self.x_shift_specific_param_name].set("0")
                    
                    self.generate_and_display(use_specific_x_shift=False, parsed_coeffs=coeffs) # Direkt mit geparsten Koeffizienten generieren
                    if self.generate_button: self.generate_button.focus_set() # Fokus setzen
                    self.update_status(f"Polynom '{poly_string}' übernommen und analysiert.")

                else:
                    messagebox.showerror("Fehler", f"Konstruktionsmethode '{target_construction_method}' nicht gefunden für Typ '{determined_type}'.")
            else:
                messagebox.showerror("Fehler", "Polynomtyp konnte nicht bestimmt werden.")

        except ValueError as e:
            messagebox.showerror("Parsing-Fehler", str(e))
            self.update_status(f"Parsing-Fehler: {str(e)[:50]}...")
        except Exception as e_gen:
            import traceback
            messagebox.showerror("Unerwarteter Fehler beim Parsen", f"{str(e_gen)}\n{traceback.format_exc()}")
            self.update_status("Unerwarteter Parsing-Fehler.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Polynomfunktionen Generator")
    parser.add_argument("--width", type=int, default=1440, help="Startbreite des Fensters")
    parser.add_argument("--height", type=int, default=970, help="Starthöhe des Fensters") # Höhe leicht erhöht
    args = parser.parse_args()

    root = tk.Tk()
    app = PolynomialAppGUI(root, initial_width=args.width, initial_height=args.height)
    root.mainloop()
