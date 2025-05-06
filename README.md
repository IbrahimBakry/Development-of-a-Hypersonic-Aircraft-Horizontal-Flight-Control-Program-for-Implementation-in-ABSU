# 🚀 **Hypersonic Aircraft Optimal Control: Fuel-Efficient Trajectory Design**  

## **📌 Project Overview**  
**Goal**: Develop an optimal control program for a **hypersonic aircraft** (Mach 4–6) during **horizontal cruise flight**, minimizing fuel consumption while meeting flight constraints. Implementable in **Automatic Onboard Control Systems (ABSU)**.  

**Key Achievements**:  
✅ **Fuel savings** via optimized trajectory planning.  
✅ **Two real-world flight paths** simulated (8,000 km & 15,000 km).  
✅ **Published 2 papers** on control algorithms and hypersonic dynamics.  

---

## **🔧 Tools & Software**  
| **Category**       | **Tools Used**                          | **Purpose** |
|--------------------|----------------------------------------|-------------|
| **Modeling**       | NASA’s GHAME model, DATCOM/SHABP       | Aerodynamic coefficients, mass properties |
| **Optimization**   | MATLAB (fmincon, NLP solver)           | Trajectory optimization, control parameter tuning |
| **Numerical Methods** | Runge-Kutta, Euler discretization   | ODE integration for dynamics |
| **Visualization**  | MATLAB Plotting, Google Earth API      | Trajectory mapping, performance graphs |

---

## **📐 Modeling & Numerical Solutions**  

### **1. 🛩️ Aircraft Dynamics Modeling**  
- **6-DOF Equations**: Derived for hypersonic flight, including:  
  - **Translational motion** (velocity, altitude, flight path angles).  
  - **Rotational dynamics** (Euler angles, angular rates).  
  - **Fuel consumption** linked to throttle control.  
- **Aerodynamic Forces**:  
  - Lift/Drag modeled as **polynomial functions** of Mach number and angle of attack (α).  
  - Side-force (Y) and moments (roll/pitch/yaw) included for **3D trajectory control**.  

### **2. ⚙️ Optimal Control Problem**  
**Objective**: Minimize fuel use → **Maximize final aircraft mass** (indirect optimization).  

#### **🔹 Pontryagin’s Maximum Principle**  
- **Hamiltonian formulation** with costate variables (Lagrange multipliers).  
- **Adjoint equations** solved numerically for optimality conditions.  
- **Challenges**: High sensitivity to initial guesses → Required **robust numerical tuning**.  

#### **🔹 Nonlinear Programming (NLP)**  
- **Interior-Point Method** (MATLAB `fmincon`) for constrained optimization.  
- **Discretization**: Trajectory split into segments → **Finite-difference approximation** of ODEs.  
- **Constraints**:  
  - Altitude = 30 km (fixed).  
  - Speed = Mach 5–6.  
  - Control limits: **α ∈ [-10°, 10°]**, **throttle ∈ [0, 1]**.  

### **3. 💻 Numerical Implementation**  
- **MATLAB Workflow**:  
  1. **Initial Guess**: Linear interpolation for states (altitude, velocity).  
  2. **Collocation Method**: Convert ODEs to algebraic constraints.  
  3. **SQP Solver**: Iteratively refine control inputs (α, μₐ, ηₜ).  
- **Convergence**: Achieved in **<50 iterations** per trajectory.  

---

## **📊 Results & Validation**  
### **Trajectory 1 (8,000 km)**  
- **Fuel saved**: ~41,027 kg (vs. non-optimal path).  
- **Control Parameters**:  
  - **α ≈ 0.7°**, **throttle ≈ 0** (nominal fuel flow).  
  - **Roll angle (μₐ) ≈ 0.5°** (minor adjustments).  

### **Trajectory 2 (15,000 km)**  
- **Fuel saved**: ~77,277 kg.  
- **Notable Behavior**:  
  - **Near-constant Mach 6** (1830 m/s).  
  - **β (sideslip) ≈ 0.05°** → Negligible impact on fuel efficiency.  

---

## **🎯 Key Takeaways**  
- **Optimal Controls**: Angle of attack (α) and throttle are **dominant** for fuel efficiency.  
- **Numerical Challenges**:  
  - **Stiff ODEs** → Adaptive step-sizing in MATLAB.  
  - **Constraint handling** → Slack variables for inequality bounds.  
- **Future Work**:  
  - Real-time **adaptive control** for turbulence.  
  - **GPU acceleration** for faster NLP convergence.  

---

## **📜 Publications**  
1. **"Fuel-Optimal Control of Hypersonic Cruise Using Pontryagin’s Principle"** (Journal of Aerospace Engineering).  
2. **"NLP-Based Trajectory Optimization for GHAME-Class Hypersonic Vehicles"** (AIAA SciTech Forum).  

---

### **🚀 Why this project?**  
This work bridges **theoretical optimal control** and **practical aerospace engineering**, enabling **long-range, fuel-efficient hypersonic travel**—critical for future aviation!  
