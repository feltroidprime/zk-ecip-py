# Import SageMath package
from sage.all import *

# Define a finite field, for example, GF(7) for field of order 7
F = GF(0x30644E72E131A029B85045B68181585D97816A916871CA8D3C208C16D87CFD47)


P.<x> = PolynomialRing(F)


# Computes (X - z)^(2^i)
def X_MIN_XP_POW_2I(x, z, i):
    return (x - z)**(2**i)

A = 0
B = 3

px = 16972146358605338978832925186803183426611046226825606566434356706141393838202
py = 13196854042874744296847103784920665732531458138601311807468784567343458201148
res = X_MIN_XP_POW_2I(x, px, 1)
print(res)

def compute_q(ri, vi, px, i):
    num = -2*ri
    den = vi
    q,r = num.quo_rem(den)
    assert r == 0, f"r = {r} != 0"
    print(f"rest: {q}")
    mod = X_MIN_XP_POW_2I(x, px, i)
    print(f"mod: {mod}")
    _, res = q.quo_rem(mod)
    print(f"res: {res}")
    
    return res

def hensel(px, py, m):
    if m>1:
        v=[]
        r=[]
        q=[]
        i=0
        v.append(py)
        print(f"v_{i}: {v[i]}")
        r_num = py**2 - (x**3 + A*x + B)
        r0, rr = r_num.quo_rem(x - px)
        assert rr == 0, f"rr = {rr} != 0"
        r.append(r0)
        print(f"r_{i}: {r[i]}")
        q.append(compute_q(r[i], v[i], px, i))
        print(f"q_{i}: {q[i]}")
        while True:
            print(f"Stage {i} check")
            check_poly = v[i]**2 - (x**3 + A*x + B)
            mod_check = X_MIN_XP_POW_2I(x, px, i)
            _, check = check_poly.quo_rem(mod_check)
            assert check == 0, f"check = {check} != 0"
            i+=1
            print(f"i: {i}")
            v.append(v[i-1] + q[i-1]*X_MIN_XP_POW_2I(x, px, i-1))

            if 2**(i-1) < m <= 2**i:
                break
        
    else:
        raise ValueError("m must be greater than 1")
    
hensel(px, py, 3)