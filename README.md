# PolyGenius

[![Python](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## 💡 Projektübersicht

**PolyGenius** ist eine interaktive Desktop-Anwendung, entwickelt mit Python und Tkinter, die primär darauf ausgelegt ist, Lehrkräften und Aufgabenstellern die **Erstellung passgenauer Polynomfunktionen für Unterrichtsmaterialien und Textaufgaben** zu erleichtern. Das Tool vereinfacht die Generierung, Analyse und Visualisierung von Polynomfunktionen bis zum 4. Grad (einschließlich biquadratischer, kubischer, quadratischer und linearer Funktionen).

Definieren Sie Funktionen durch spezifische Eigenschaften wie Nullstellen, Extrema, Wendepunkte, Skalierungsfaktoren oder direkt über ihre Koeffizienten. PolyGenius berechnet das entsprechende Polynom, zeigt seine Ableitungen sowie charakteristische Punkte an und visualisiert die Funktion in einem interaktiven Graphen.

**Neu:** Geben Sie Polynome direkt als Text ein (z.B. "x^3 - 2x + 4")! PolyGenius erkennt den Typ und die Koeffizienten automatisch und stellt die Funktion zur weiteren Analyse bereit.

Mit Optionen wie dem **Erzwingen ganzzahliger Koeffizienten**, dem **"Ich fühle mich glücklich!"-Modus** für zufällige, aufgabenfreundliche Polynome und der Möglichkeit, eine **detaillierte Kurvendiskussion** (inklusive Symmetrie, Monotonie und Krümmung) zu exportieren, ist PolyGenius ein vielseitiges Werkzeug. Es unterstützt nicht nur bei der Erstellung von Lehrmaterial, sondern dient auch Schülern und Studenten als Lernhilfe zum besseren Verständnis von Polynomfunktionen.

---

## ✨ Features

* **Vielseitige Polynomtypen:** Erstellen Sie Funktionen bis zum 4. Grad:
    * **Biquadratische Polynome** ($Ax^4 + Bx^2 + D$)
    * **Kubische Polynome** ($ax^3 + bx^2 + cx + d$)
    * **Quadratische Polynome** ($bx^2 + cx + d$)
    * **Lineare Polynome** ($cx + d$)
* **Flexible Konstruktionsmethoden:**
    * **Direkte Eingabe:** Geben Sie die Polynomfunktion als String ein (z.B. "2x^3 - x^2 + 5x - 1") oder definieren Sie sie über ihre Koeffizienten.
    * **Kubische Polynome:** Definieren Sie über Sattelpunkte, Extrema-Lagen (x/y-Werte), Wendepunkt & Abstand zu Extrema, oder Nullstellen.
    * **Biquadratische Polynome:** Konstruktion aus äußeren und zentralen Extrema, Wendepunkten und zentralem Extremum, oder direkte Koeffizienteneingabe ($A, B, D$).
    * **Quadratische Polynome:** Konstruktion aus Scheitelpunkt und Skalierungsfaktor oder direkte Koeffizienteneingabe.
    * **Lineare Polynome:** Aus Nullstelle und y-Achsenabschnitt oder direkte Koeffizienteneingabe.
    * **X-Verschiebung:** Optionale Verschiebung des Basispolynoms entlang der x-Achse für kubische und biquadratische Funktionen.
* **Präzise Berechnungen:** Verwendet `fractions` für exakte Koeffizienten und `numpy` für robuste Nullstellensuche.
* **Detaillierte Analyse & Export:**
    * Zeigt das (ggf. verschobene) Polynom $P(x)$, seine erste $P'(x)$, zweite $P''(x)$ und dritte $P'''(x)$ Ableitung an.
    * Berechnet Nullstellen, kritische Stellen (Extrema) mit Klassifizierung (Hoch-, Tief-, Sattelpunkt), und Wendestellen.
    * **NEU:** Exportieren Sie eine vollständige Kurvendiskussion als Text, inklusive Symmetrieanalyse, Monotonie- und Krümmungsintervallen.
* **Interaktiver Graph:**
    * Visualisiert das generierte Polynom zusammen mit seinen Nullstellen, Extrema und Wendepunkten.
    * Matplotlib-Toolbar für Zoom, Verschiebung und Speichern des Graphen.
    * Optimierte automatische Achsenlimits mit der Möglichkeit, manuell herauszuzoomen, um den Graphen über ein weites Intervall zu betrachten.
* **Benutzerfreundliche GUI:**
    * Intuitive Oberfläche mit Dropdown-Menüs und direkter Eingabe.
    * Statusleiste für Rückmeldungen.
    * "Ich fühle mich glücklich!"-Modus für zufällige Polynome.
    * Option zum Erzwingen ganzzahliger Koeffizienten.
    * Startfenstergröße über Kommandozeilenargumente anpassbar.

---

### Voraussetzungen

* Python 3.x
* Die folgenden Python-Bibliotheken:
    * `numpy`
    * `matplotlib`

---

## 🚀 Erste Schritte

1.  **Klonen oder Herunterladen:**
    ```bash
    git clone [https://github.com/FisherSteve/PolyGenius.git](https://github.com/FisherSteve/PolyGenius.git) 
    cd PolyGenius
    ```
    

2.  **Abhängigkeiten installieren (falls noch nicht geschehen):**
    ```bash
    pip install numpy matplotlib
    ```

3.  **Anwendung starten:**
    Führe das Skript über die Kommandozeile aus:
    ```bash
    python PolyGeniusBeta.py
    ```
    
    Optional kannst du die Startfenstergröße übergeben:
    ```bash
    python dein_programm_name.py --width 1600 --height 1000 --scale 0.6
    ```

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
    Navigiere in der Kommandozeile zu dem Ordner, in dem sich deine Python-Datei (z.B. `poly_genius.py`) befindet. Zum Beispiel:
    ```bash
    cd C:\Pfad\zu\deinem\Projektordner
    ```

3.  **Die .exe-Datei erstellen:**
    Führe den folgenden Befehl in der Kommandozeile aus. Er teilt PyInstaller mit, dass eine einzelne, fensterbasierte (`--windowed`) `.exe`-Datei erstellt werden soll:
    ```bash
    pyinstaller --onefile --windowed PolyGeniusBeta.py
    ```
    * `--onefile`: Packt die gesamte Anwendung in **eine einzige `.exe`-Datei**.
    * `--windowed` (oder `--noconsole`): Verhindert, dass beim Start der `.exe` ein zusätzliches Konsolenfenster geöffnet wird.
    * *(Optional: Für ein eigenes Icon füge `--icon=dein_icon.ico` hinzu.)*

4.  **Die fertige .exe finden:**
    Nachdem PyInstaller seine Arbeit beendet hat, findest du die generierte `.exe`-Datei im Ordner `dist` in deinem Projektverzeichnis.

---

## Lizenz

Dieses Projekt ist unter der MIT-Lizenz lizenziert.

