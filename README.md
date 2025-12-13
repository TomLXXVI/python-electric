# python-electric

**python-electric** is a lightweight Python package designed to support the conceptual design and verification of 
low-voltage (LV) electrical installations *without requiring manufacturer-specific time–current curves*. The library 
implements the fundamental protection principles from **IEC 60364**, **IEC 60898-1**, and **IEC 60947-2**, enabling 
engineers to perform:

* Cable sizing
* Overload protection checks
* Short-circuit protection checks
* Basic current-based selectivity analysis
* Time-based selectivity estimation for industrial circuit breakers
* Earthing system verification (TN) including automatic disconnection of supply (ADS)

This README provides an overview of the package, its design philosophy, main features, limitations, and typical 
workflows.

---

## 1. Design Philosophy

The goal of *python-electric* is not to replace full professional software (e.g., Ecodial, Caneco, Simaris), but to 
provide:

* A **transparent, normative, calculation-driven** approach
* A tool that works **without manufacturer curves**
* A minimal but robust model of circuit breakers, cables, and LV networks
* A framework suitable for education, prototyping, and conceptual design

Where detailed manufacturer data is unavailable or unnecessary, the package uses **normative parameters**:

* Conventional non-tripping current ($I_{nf}$)
* Conventional tripping current ($I_{f}$)
* Magnetic trip current band ($I_{m,min}$ – $I_{m,max}$)
* Maximum magnetic tripping time ($t_{m,lim}$)
* Ultimate breaking capacity ($I_{cu}$)
* Let-through energy ($I^2t$)

This allows the user to verify **overload protection**, **thermal short‑circuit withstand**, and **basic selectivity**, 
in accordance with IEC principles.

---

## 2. Core Components

### **2.1 Circuit Breakers**

Circuit breakers are modeled with their essential IEC parameters. The package distinguishes between:

* **Residential MCBs** (IEC 60898-1)

  * Fixed time–current behaviour
  * No adjustable delay in the magnetic region
  * Maximum magnetic clearing time typically ≤ 100 ms

* **Industrial MCCBs/ACBs** (IEC 60947-2)

  * Optional adjustable short-time delay
  * Support for time-based selectivity

The class exposes flags such as:

* `standard` (RESIDENTIAL vs INDUSTRIAL)
* `category` (B, C, D, ADJUSTABLE)
* `has_adjustable_delay` (derived property)

---

### **2.2 Cables**

Cables are characterized by:

* Cross-sectional area (S)
* Material and insulation type (affecting the *k*-factor)
* Ampacity (I_z)
* Joule integral (k²·S²)
* Short-circuit currents at upstream and downstream ends

Cables provide methods to evaluate:

* Thermal withstand under short-circuit
* Overload capacity
* TN earthing ADS compliance (fault loop impedance + touch voltage)

---

## 3. Protection Checks

### **3.1 Overload Protection**

The package verifies the IEC criteria:

$$
I_b \le I_n \le I_z
$$

$$
I_2 \le 1.45 I_z
$$

where I₂ is the conventional tripping current.

---

### **3.2 Short-Circuit Protection**

Short-circuit verification includes:

* **Breaking capacity check:**
  $$
  I_{k,\text{max}} \le I_{cu}
  $$

* **Thermal protection of the cable:**
  $$
  I^2 t_{\text{device}} \le k^2 S^2
  $$

* **Magnetic pickup verification** for minimum short-circuit currents.

---

## 4. Selectivity

### **4.1 Current-Based Selectivity**

The upstream and downstream circuit breakers are compared in terms of:

* Their thermal characteristics (I₁/I₂)
* Magnetic trip bands (I_m,min / I_m,max)

If the characteristic bands do not overlap, **full current-based selectivity** is achieved.

---

### **4.2 Time-Based Selectivity for Industrial Breakers**

For **industrial, adjustable** breakers (IEC 60947-2), the package estimates the **maximum allowable delay** for the 
upstream breaker, based on:

* Thermal withstand of the upstream cable
* Maximum allowed fault duration from ADS (TN, BB2 default)
* Short-circuit current levels in both circuits

Output:

* `True` → full selectivity
* `False` → no selectivity possible
* `Quantity (time)` → maximum permissible time delay (upper bound)

**Note:** Residential MCBs do not support adjustable delay; for them, no time-based value is returned.

---

## 5. Earthing and Fault Protection

Currently, the package supports **TN systems**, including:

* Automatic disconnection of supply (ADS)
* Fault current calculation
* Touch-voltage check against IEC safety curves (BB2 by default)

Future versions may introduce TT or IT support.

---

## 6. Example Workflow

An example Jupyter notebook is provided (`ex4_installation_sizing_01.ipynb`) demonstrating:

1. Defining transformers, cables, and loads
2. Selecting circuit breakers
3. Computing short-circuit currents
4. Checking overload & short-circuit protection
5. Checking ADS compliance
6. Evaluating selectivity

This notebook reproduces a worked example from the reference book *Laagspanningsinstallaties: technologie en ontwerp*.

---

## 7. Limitations

* The model intentionally avoids the use of manufacturer-specific curves.
* Selectivity is computed only in a **normative, conservative** sense.
* Residential MCBs cannot perform time-based selectivity.
* Arc-energy reduction, breaker coordination classes, and cascading effects (back-up protection) are outside scope.
* The package currently supports only **TN earthing systems**.

---

## 8. Roadmap

Planned extensions include:

* TT/IT earthing support
* Improved selectivity engine with manufacturer curve import (optional)
* Visualization utilities for t–I envelopes
* Report generation for compliance summaries

---

## 9. License

To be determined.

---

## 10. Acknowledgements

This package was developed as a practical engineering tool inspired by the IEC 60364 family and the textbook 
*Laagspanningsinstallaties: technologie en ontwerp*. It aims to give engineers and students a transparent, 
calculation-first approach to LV design.

---

*This is a preliminary version of the README. It will evolve as the package grows.*
