# Import SageMath package
from sage.all import *

# bn254
F = GF(0x30644E72E131A029B85045B68181585D97816A916871CA8D3C208C16D87CFD47)


P.<x> = PolynomialRing(F)



A = F(0)
B = F(3)

px = F(16972146358605338978832925186803183426611046226825606566434356706141393838202)
py = F(13196854042874744296847103784920665732531458138601311807468784567343458201148)

def is_on_curve(px, py):
    if py**2 == px**3 + A*px + B:
        print(f"Point ({px}, {py}) is on the curve.")
        return True
    else:
        print(f"Point ({px}, {py}) is not on the curve.")
        return False


assert is_on_curve(px, py)

# Computes (X - z)^(2^i)
def X_MIN_XP_POW_2I(x, z, i):
    return (x - z)**(2**i)

def compute_q(ri, vi, px, i):
    num = -2*ri
    den = vi
    q,r = num.quo_rem(den)
    assert r == 0, f"r = {r} != 0"
    print(f"rest: {q}")
    mod = X_MIN_XP_POW_2I(x, px, i)
    _, res = q.quo_rem(mod)
    print(f"res: {res}")
    
    return res

def check_stage(i, vi, xp):
    check_poly = vi**2 - (x**3 + A*x + B)
    mod_check = X_MIN_XP_POW_2I(x, xp, i)
    _, check = check_poly.quo_rem(mod_check)
    assert check == 0, f"check = {check} != 0"
    print(f"Stage {i} working, (x-xp)^(2^{i}) does divide v_{i}^2 - (x^3 + A*x + B)")
    return True

def hensel(xp, yp, m):
    if m>1:
        v=[]
        r=[]
        q=[]
        i=0
        v.append(yp)
        r_num = yp**2 - (x**3 + A*x + B)
        r0, rr = r_num.quo_rem(x - xp)
        assert rr == 0, f"rr = {rr} != 0"
        r.append(r0)
        print(f"r_{i}: {r[i]}")
        q.append(compute_q(r[i], v[i], xp, i))
        print(f"q_{i}: {q[i]}")

        check_stage(i, v[i], xp)

        while True:
            i+=1
            print(f"i: {i}")
            v.append(v[i-1] + q[i-1]*X_MIN_XP_POW_2I(x, px, i-1))

            check_stage(i, v[i], px)
            # Todo: compute r_i+1, q_i+1 after check stage works for i == 1 
            if 2**(i-1) < m <= 2**i:
                break
        
    else:
        raise ValueError("m must be greater than 1")
    
hensel(px, py, 3)