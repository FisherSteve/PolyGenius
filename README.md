# PolyGenius

[![Python](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## 💡 Projektübersicht

**PolyGenius** ist eine interaktive Desktop-Anwendung, die mit Python und Tkinter entwickelt wurde, um die Generierung, Analyse und Visualisierung von Polynomfunktionen (kubisch, quadratisch, linear) zu vereinfachen. Egal, ob Sie bestimmte Eigenschaften wie Nullstellen, Extrema, Wendepunkte oder Skalierungsfaktoren vorgeben möchten – PolyGenius berechnet das entsprechende Polynom, zeigt seine Ableitungen und kritischen Punkte an und visualisiert es in einem interaktiven Graphen.

Mit der Option, **ganzzahlige Koeffizienten zu erzwingen** und dem **"Ich fühle mich glücklich!"-Modus** für zufällige Polynome, ist PolyGenius ein ideales Werkzeug für Schüler, Studenten und alle, die ein tieferes Verständnis von Polynomfunktionen entwickeln möchten.

---

## ✨ Features

* **Vielseitige Polynomtypen:** Erstellen Sie kubische, quadratische und lineare Funktionen.
* **Flexible Konstruktionsmethoden:**
    * **Kubische Polynome:** Definieren Sie über Sattelpunkte, Extrema-Lagen (x/y-Werte), Wendepunkt & Abstand zu Extrema, oder Nullstellen.
    * **Quadratische Polynome:** Konstruktion aus Scheitelpunkt und Skalierungsfaktor oder direkte Koeffizienteneingabe.
    * **Lineare Polynome:** Aus Nullstelle und y-Achsenabschnitt oder direkte Koeffizienteneingabe.
* **Präzise Berechnungen:** Verwendet `fractions` für exakte Koeffizienten und `numpy` für robuste Nullstellensuche.
* **Detaillierte Analyse:** Zeigt das Polynom P(x), seine erste P'(x) und zweite P''(x) Ableitung, Nullstellen, kritische Stellen (Extrema), Wendestellen und Sattelpunktinformationen an.
* **Interaktiver Graph:** Visualisiert das generierte Polynom zusammen mit seinen Nullstellen, Extrema und Wendepunkten.
* **"Ich fühle mich glücklich!"-Modus:** Generiert zufällige Polynome und deren Analyse auf Knopfdruck.
* **Anpassbare Koeffizienten:** Option zum Erzwingen ganzzahliger Koeffizienten für bestimmte Konstruktionsmethoden.
* **Benutzerfreundliche GUI:** Intuitive Oberfläche, die auch für Nicht-Programmierer leicht zu bedienen ist.

---

## 🚀 Erste Schritte

### Voraussetzungen

* Python 3.x
* Die folgenden Python-Bibliotheken:
    * `tkinter` (standardmäßig in Python enthalten)
    * `numpy`
    * `matplotlib`

Installieren Sie die benötigten Bibliotheken mit pip:

```bash
pip install numpy matplotlib
