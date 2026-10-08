"""Exact coefficient certificate for higher-crossing all-spin remainders.

This proves only the stated algebraic crossing-kernel sign, not a dual bound.
All stored floating multipliers are interpreted as their exact binary rationals.
"""
from __future__ import annotations
import argparse, hashlib, json, math, sys
from fractions import Fraction
from pathlib import Path
import numpy as np
import sympy as sy

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'research'))
from crossing_kernels import crossing_catalogue, build_crossing_bank


def coefficient_polynomial(witness):
    data = np.load(witness)
    multipliers = dict(zip(map(str, data['names']), data['multipliers']))
    degree = max(int(name[1:].split('_')[0]) for name in multipliers if name.startswith('X'))
    k, z = sy.symbols('k z')
    c = sy.Rational(4, 3)
    expr = sy.S.Zero
    for row in crossing_catalogue(degree):
        if row['name'] not in multipliers:
            continue
        lam = sy.Rational(float(multipliers[row['name']]))
        for j, beta in enumerate(row['beta']):
            falling = sy.prod(k - i for i in range(j)) / sy.factorial(j)
            expr += lam * beta * falling / c**j * z**(row['degree'] - 2 - j)
    return sy.Poly(sy.expand(expr), k, z)


def bernstein_univariate(expr, variable, a, b, degree=None):
    t = sy.Symbol('t')
    poly = sy.Poly(sy.expand(expr.subs(variable, a + (b-a)*t)), t)
    minimum_degree = max(0, poly.degree()) if not poly.is_zero else 0
    degree = minimum_degree if degree is None else degree
    if degree < minimum_degree:
        raise ValueError('Insufficient Bernstein degree')
    return [sy.factor(sum(poly.nth(j)*sy.binomial(i,j)/sy.binomial(degree,j)
                         for j in range(i+1))) for i in range(degree+1)]


def optimizer_guard_matrix(names, cutoff=10, degree=8):
    """Linear exact guards A@lambda<=bounds proving F(k,z)<=-eta for k>=K.

    Returns rational row entries and the power_y for each row. Choose bounds
    -eta on power_y==0, zero on all remaining rows. These are sufficient
    search restrictions, not extra assumptions on physical amplitudes.
    """
    names=list(map(str,names)); k,z,y=sy.symbols('k z y');c=sy.Rational(4,3)
    by_name={row['name']:row for row in crossing_catalogue(degree)}
    n_y=degree-2;n_z=degree-3
    polynomials=[]
    for name in names:
        if name not in by_name:
            polynomials.append(sy.Poly(0,y));continue
        row=by_name[name]
        expr=sum(beta*sy.prod(k-i for i in range(j))/sy.factorial(j)/c**j*z**(row['degree']-2-j)
                 for j,beta in enumerate(row['beta']) if beta)
        polynomials.append(sy.Poly(sy.expand(expr.subs(k,cutoff+y)),y))
    matrix=[];powers=[]
    for power in range(n_y+1):
        columns=[bernstein_univariate(poly.nth(power),z,sy.S.Zero,sy.Rational(3,8),n_z)
                 for poly in polynomials]
        for index in range(n_z+1):
            row=[column[index] for column in columns]
            if any(row):
                matrix.append(row);powers.append(power)
    return matrix,powers


def shifted_certificate(poly, cutoff, subdivisions=1):
    k, z = poly.gens
    y = sy.Symbol('y')
    shifted = sy.Poly(sy.expand(poly.as_expr().subs(k, cutoff+y)), y)
    records = []
    for i in range(shifted.degree()+1):
        for interval in range(subdivisions):
            a = sy.Rational(3*interval, 8*subdivisions)
            b = sy.Rational(3*(interval+1), 8*subdivisions)
            coefficients = bernstein_univariate(shifted.nth(i), z, a, b)
            if max(coefficients) > 0:
                return None
            records.append({'power_y': i, 'z_interval': [str(a), str(b)],
                            'bernstein_coefficients': list(map(str, coefficients))})
    return records


def positive_low_coefficients(poly, cutoff, subdivisions=16):
    k, z = poly.gens
    out = []
    for i in range(cutoff):
        bounds = []
        for interval in range(subdivisions):
            a = sy.Rational(3*interval, 8*subdivisions)
            b = sy.Rational(3*(interval+1), 8*subdivisions)
            coeff = bernstein_univariate(poly.as_expr().subs(k, i), z, a, b)
            bounds.append(max(coeff))
        upper = max(bounds)
        out.append({'k': i, 'global_upper_rational': str(upper),
                    'global_upper_approx': float(upper),
                    'nonpositive_for_all_energies': bool(upper <= 0)})
    return out


def replay_identity(poly, witness):
    """Independent mpmath Legendre-series replay against the existing Jacobi code."""
    import mpmath as mp
    mp.mp.dps=80
    data=np.load(witness)
    multipliers=dict(zip(map(str,data['names']),data['multipliers']))
    k,z=poly.gens
    energies=[4.01,4.5,8.,20.,100.,1000.,10000.]
    spins=[0,2,4,8,10,20,40,80]
    degree=max(int(name[1:].split('_')[0]) for name in multipliers if name.startswith('X'))
    bank=build_crossing_bank(energies,spins,degree,normalized=True)
    reference=sum(float(multipliers.get(name,0))*row for name,row in zip(bank['labels'],bank['v']))
    worst=0.
    for i,ell in enumerate(spins):
        for j,energy in enumerate(energies):
            s=mp.mpf(energy);delta=s-4;Q=s-mp.mpf(4)/3;x=mp.mpf(4)/(3*delta)
            coefficients=[mp.mpf(1)]
            for n in range(1,ell+1):
                coefficients.append(coefficients[-1]*(ell*(ell+1)-(n-1)*n)/n**2*x)
            numerator=mp.mpf(0)
            for n,b in enumerate(coefficients):
                expression=poly.as_expr().subs(k,n)
                rational_terms=sy.Poly(expression,z)
                value=mp.mpf(0)
                for (power,),coefficient in rational_terms.terms():
                    value+=mp.mpf(int(coefficient.p))/int(coefficient.q)/Q**power
                numerator+=b*value
            actual=float(numerator/sum(coefficients))
            relative=abs(actual-reference[i,j])/(1+abs(actual))
            worst=max(worst,relative)
    if worst>2e-10:
        raise RuntimeError('Crossing-series identity replay mismatch '+str(worst))
    return {'arithmetic':'80-digit mpmath positive Legendre expansion vs existing scipy Jacobi kernels',
            'energies':energies,'spins':spins,'max_relative_scale_error':worst,
            'scope':'Diagnostic implementation cross-check; the sign certificate itself uses exact rational arithmetic'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('witness', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    poly = coefficient_polynomial(args.witness)
    certificate = None
    for cutoff in range(2, 65):
        for subdivisions in [1, 2, 4, 8, 16]:
            certificate = shifted_certificate(poly, cutoff, subdivisions)
            if certificate is not None:
                break
        if certificate is not None:
            break
    if certificate is None:
        raise RuntimeError('Sufficient shifted-power/Bernstein sign test failed')
    negative_margin=-max(sy.Rational(value) for record in certificate if record['power_y']==0
                         for value in record['bernstein_coefficients'])
    result = {
        'scope': 'Exact all-energy crossing-series sign certificate; NOT a full dual bound',
        'witness': str(args.witness.resolve()),
        'witness_sha256': hashlib.sha256(args.witness.read_bytes()).hexdigest(),
        'identity': 'Xcombined/p=sum(B_k*x^k*F(k,1/Q))/sum(B_k*x^k), x=4/(3*(mu-4))',
        'polynomial': str(poly.as_expr()),
        'cutoff_k': cutoff,
        'energy_domain': '4<mu<infinity; z=1/(mu-4/3) in (0,3/8)',
        'conclusion': 'F(k,z)<=0 for every real k>=cutoff_k and every z in [0,3/8]',
        'uniform_negative_margin_rational': str(negative_margin),
        'uniform_negative_margin_approx': float(negative_margin),
        'stronger_conclusion': 'F(k,z)<=-uniform_negative_margin_rational on the same domain',
        'certificate': certificate,
        'low_coefficients': positive_low_coefficients(poly, cutoff),
        'identity_replay': replay_identity(poly,args.witness),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2), encoding='utf8')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('certificate','polynomial')}, indent=2))


if __name__ == '__main__':
    main()
