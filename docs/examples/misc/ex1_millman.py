import math

from python_electric import Q_
from python_electric import calc

U_eff_line = Q_(400, 'V').m
P_l1 = Q_(400, 'W').m
P_l2 = Q_(100, 'W').m
P_l3 = Q_(1e-12, 'W').m


U_eff_ph = U_eff_line / math.sqrt(3)

I_eff_1 = P_l1 / U_eff_ph
I_eff_2 = P_l2 / U_eff_ph
I_eff_3 = P_l3 / U_eff_ph

I_1 = calc.phasor(I_eff_1, 0.0)
I_2 = calc.phasor(I_eff_2, -120.0)
I_3 = calc.phasor(I_eff_3, -240.0)

U_1 = calc.phasor(U_eff_ph, 0.0)
U_2 = calc.phasor(U_eff_ph, -120.0)
U_3 = calc.phasor(U_eff_ph, -240.0)

R_l1 = U_1 / I_1
R_l2 = U_2 / I_2
R_l3 = U_3 / I_3


millman = calc.MillmanTheorem(
    phase_voltages=(U_1, U_2, U_3),
    source_impedances=None,
    load_impedances=(R_l1, R_l2, R_l3),
    neutral_impedance=float("inf")
)

print(f"deltaU_neutral = {calc.polar(millman.dU_neutral)}")

I_line = millman.I_line
print(f"I_1 = {calc.polar(I_line[0])}")
print(f"I_2 = {calc.polar(I_line[1])}")
print(f"I_3 = {calc.polar(I_line[2])}")

U_l1 = I_line[0] * R_l1
U_l2 = I_line[1] * R_l2
U_l3 = I_line[2] * R_l3
print(f"U_l1 = {calc.polar(U_l1)}")
print(f"U_l2 = {calc.polar(U_l2)}")
print(f"U_l3 = {calc.polar(U_l3)}")
