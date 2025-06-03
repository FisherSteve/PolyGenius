# PolyGenius

[![Python](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## 💡 Projektübersicht

**PolyGenius** ist eine interaktive Desktop-Anwendung, die mit Python und Tkinter entwickelt wurde, um die Generierung, Analyse und Visualisierung von Polynomfunktionen (kubisch, quadratisch, linear) zu vereinfachen. Egal, ob Sie bestimmte Eigenschaften wie Nullstellen, Extrema, Wendepunkte oder Skalierungsfaktoren vorgeben möchten – PolyGenius berechnet das entsprechende Polynom, zeigt seine Ableitungen und kritischen Punkte an und visualisiert es in einem interaktiven Graphen.

Mit der Option, **ganzzahlige Koeffizienten zu erzwingen** und dem **"Ich fühle mich glücklich!"-Modus** für zufällige Polynome, ist PolyGenius ein ideales Werkzeug für Schüler, Studenten und alle, die ein tieferes Verständnis von Polynomfunktionen entwickeln möchten oder Lehrer, welche Aufgaben erstellen wollen.

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


### Voraussetzungen

* Python 3.x
* Die folgenden Python-Bibliotheken:
    * `tkinter` (standardmäßig in Python enthalten)
    * `numpy`
    * `matplotlib`


---

## 📦 Als ausführbare Datei (.exe) verteilen

Möchtest du **PolyGenius** an Freunde oder Kollegen weitergeben, die keine Python-Umgebung auf ihrem Computer haben? Kein Problem! Du kannst ganz einfach eine eigenständige `.exe`-Datei erstellen, die auf Windows-Systemen direkt ausführbar ist. Dafür nutzen wir **PyInstaller**, ein beliebtes Tool, das dein Python-Skript und alle Abhängigkeiten in eine einzige ausführbare Datei verpackt.

### Schritt-für-Schritt-Anleitung

1.  **PyInstaller installieren:**
    Öffne deine Kommandozeile (CMD oder PowerShell) und installiere PyInstaller, falls du es noch nicht hast:
    ```bash
    pip install pyinstaller
    ```

2.  **In das Projektverzeichnis wechseln:**
    Navigiere in der Kommandozeile zu dem Ordner, in dem sich deine Datei `PolyGeniusBeta.py` befindet. Zum Beispiel:
    ```bash
    cd C:\Users\DeinNutzername\Dokumente\DeinProjektordner
    ```
    *(Passe den Pfad entsprechend an!)*

3.  **Die .exe-Datei erstellen:**
    Führe den folgenden Befehl in der Kommandozeile aus. Er teilt PyInstaller mit, dass eine einzelne, fensterbasierte (`--windowed`) `.exe`-Datei aus deiner `PolyGeniusBeta.py` erstellt werden soll:
    ```bash
    pyinstaller --onefile --windowed PolyGeniusBeta.py
    ```
    * `--onefile`: Packt die gesamte Anwendung in **eine einzige `.exe`-Datei**. Das ist super praktisch für die Verteilung.
    * `--windowed` (oder `--noconsole`): Verhindert, dass beim Start der `.exe` ein zusätzliches schwarzes Konsolenfenster im Hintergrund geöffnet wird – ideal für Tkinter-GUIs.
    * *(Optional: Möchtest du ein eigenes Icon für deine `.exe`? Füge `--icon=mein_icon.ico` hinzu, wobei `mein_icon.ico` der Pfad zu deiner Icon-Datei im `.ico`-Format ist.)*

4.  **Die fertige .exe finden:**
    Nachdem PyInstaller seine Arbeit beendet hat (das kann ein paar Minuten dauern), findest du die generierte `.exe`-Datei im Ordner `dist` in deinem Projektverzeichnis.
    
    Du kannst sie dann an deine Kollegen weitergeben!

---
