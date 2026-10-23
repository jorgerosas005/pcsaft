import pcsaft

c = pcsaft.get("ethane")
print(pcsaft.spinodal(250.0, c))   # (eta_v, eta_l, Pmax, Pmin))
print(pcsaft.spinodal(300.0, c))   # the loop has shrunk: ~39 to ~45 bar
print(pcsaft.spinodal(320.0, c))   # None: above Tc, no loop at all