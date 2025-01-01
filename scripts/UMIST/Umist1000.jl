### Packages ###
using Catalyst
using DifferentialEquations
using Plots
using Sundials, LSODA
using OrdinaryDiffEq
using Symbolics
using DiffEqDevTools
using ODEInterface, ODEInterfaceDiffEq
using ModelingToolkit


### Timespan ###
seconds_per_year = 3600 * 24 * 365
tspan = (0.0, 30000 * seconds_per_year) # ~30 thousand yrs


### Initial Conditions ###
u0 = [
      0.0,   # 1: C⁻
      0.0,   # 2: C
      0.0,   # 3: C2
      0.0,   # 4: e
      0.0,   # 5: CH2
      0.0,   # 6: C2H2
      0.0,   # 7: CH
      0.0,   # 8: C2H
      0.0,   # 9: CO2
      0.0,   # 10: CO
      0.0,   # 11: H2O
      0.0,   # 12: H2CO
      0.0,   # 13: N
      0.0,   # 14: CN
      0.0,   # 15: NH
      0.0,   # 16: HCN
      0.0,   # 17: O2
      0.0,   # 18: O
      0.0,   # 19: OH
      0.0,   # 20: HCO
      0.0,   # 21: C2⁻
      0.0,   # 22: C4
      0.0,   # 23: C3
      0.0,   # 24: C5
      0.0,   # 25: C6
      0.0,   # 26: C7
      0.0,   # 27: C8
      0.0,   # 28: C9
      0.0,   # 29: C10
      0.0,   # 30: C3⁻
      0.0,   # 31: C4⁻
      0.0,   # 32: C5⁻
      0.0,   # 33: C6⁻
      0.0,   # 34: C7⁻
      0.0,   # 35: C8⁻
      0.0,   # 36: C10⁻
      0.0,   # 37: C2H⁻
      0.0,   # 38: C3H
      0.0,   # 39: C3H⁻
      0.0,   # 40: C4H
      0.0,   # 41: C3N⁻
      0.0,   # 42: C4N
      0.0,   # 43: C4H⁻
      0.0,   # 44: C5H
      0.0,   # 45: C5H⁻
      0.0,   # 46: C6H
      0.0,   # 47: C6H⁻
      0.0,   # 48: C7H
      0.0,   # 49: C7H⁻
      0.0,   # 50: C8H
      0.0,   # 51: C8H⁻
      0.0,   # 52: C9H
      0.0,   # 53: C9⁻
      0.0,   # 54: C9H⁻
      0.0,   # 55: C10H
      0.0,   # 56: CH⁻
      0.0,   # 57: O⁻
      0.0,   # 58: OH⁻
      0.0,   # 59: S⁻
      0.0,   # 60: CS
      0.0,   # 61: CH3
      0.0,   # 62: CN⁻
      0.0,   # 63: CH3CN
      0.0,   # 64: CH3OH
      0.0,   # 65: HCO⁺
      0.0,   # 66: OCS
      0.0,   # 67: H⁻
      0.0,   # 68: CH4
      0.0,   # 69: H
      0.0,   # 70: H2
      0.0,   # 71: NH2
      0.0,   # 72: NH3
      0.0,   # 73: C10H⁻
      0.0,   # 74: C10H2
      0.0,   # 75: H2CCC
      0.0,   # 76: HC3N
      0.0,   # 77: HC4H
      0.0,   # 78: C5H2
      0.0,   # 79: C5N⁻
      0.0,   # 80: HC5N
      0.0,   # 81: C6H2
      0.0,   # 82: C7H2
      0.0,   # 83: C8H2
      0.0,   # 84: C9H2
      0.0,   # 85: HS
      0.0,   # 86: C2N
      0.0,   # 87: C3N
      0.0,   # 88: C5N
      0.0,   # 89: C7N
      0.0,   # 90: HC7N
      0.0,   # 91: C9N
      0.0,   # 92: HC9N
      0.0,   # 93: NO
      0.0,   # 94: NS
      0.0,   # 95: NO2
      0.0,   # 96: SO2
      0.0,   # 97: C2O
      0.0,   # 98: SO
      0.0,   # 99: H⁺
      0.0,   # 100: HNC
      0.0,   # 101: HOC⁺
      0.0,   # 102: O2⁻
      0.0,   # 103: C⁺
      0.0,   # 104: C10H⁺
      0.0,   # 105: C2H4
      0.0,   # 106: C2H4⁺
      0.0,   # 107: C2H5
      0.0,   # 108: C2H5⁺
      0.0,   # 109: C2H5OH
      0.0,   # 110: C2H5OH⁺
      0.0,   # 111: C2O⁺
      0.0,   # 112: C2S
      0.0,   # 113: C2S⁺
      0.0,   # 114: C3O
      0.0,   # 115: C3O⁺
      0.0,   # 116: C3S
      0.0,   # 117: C3S⁺
      0.0,   # 118: C4H3
      0.0,   # 119: C4H3⁺
      0.0,   # 120: C4S
      0.0,   # 121: C4S⁺
      0.0,   # 122: C6H6
      0.0,   # 123: C6H6⁺
      0.0,   # 124: CCP
      0.0,   # 125: CCP⁺
      0.0,   # 126: CCl
      0.0,   # 127: CCl⁺
      0.0,   # 128: CH2⁺
      0.0,   # 129: CH2CCH2
      0.0,   # 130: C3H4⁺
      0.0,   # 131: CH2CCH
      0.0,   # 132: CH2CCH⁺
      0.0,   # 133: CH2CN
      0.0,   # 134: CH2CN⁺
      0.0,   # 135: CH2CO
      0.0,   # 136: CH2CO⁺
      0.0,   # 137: CH3CCH
      0.0,   # 138: CH3CHCH2
      0.0,   # 139: C3H6⁺
      0.0,   # 140: CH3CHO
      0.0,   # 141: CH3CHO⁺
      0.0,   # 142: CH3COCH3
      0.0,   # 143: CH3COCH3⁺
      0.0,   # 144: CH3OCH3
      0.0,   # 145: CH3OCH3⁺
      0.0,   # 146: CH⁺
      0.0,   # 147: CP
      0.0,   # 148: CP⁺
      0.0,   # 149: ClO
      0.0,   # 150: ClO⁺
      0.0,   # 151: Fe
      0.0,   # 152: Fe⁺
      0.0,   # 153: H2CO⁺
      0.0,   # 154: H2S
      0.0,   # 155: H2S⁺
      0.0,   # 156: H2SiO
      0.0,   # 157: H2SiO⁺
      0.0,   # 158: C4H2⁺
      0.0,   # 159: HCOOCH3
      0.0,   # 160: COOCH4⁺
      0.0,   # 161: HCP
      0.0,   # 162: HCP⁺
      0.0,   # 163: HPO
      0.0,   # 164: HPO⁺
      0.0,   # 165: Mg
      0.0,   # 166: Mg⁺
      0.0,   # 167: NCCN
      0.0,   # 168: C2N⁺
      0.0,   # 169: CNC⁺
      0.0,   # 170: NH3⁺
      0.0,   # 171: NO⁺
      0.0,   # 172: NS⁺
      0.0,   # 173: Na
      0.0,   # 174: Na⁺
      0.0,   # 175: OCS⁺
      0.0,   # 176: P
      0.0,   # 177: P⁺
      0.0,   # 178: PH
      0.0,   # 179: PH⁺
      0.0,   # 180: PO
      0.0,   # 181: PO⁺
      0.0,   # 182: SO⁺
      0.0,   # 183: Si
      0.0,   # 184: Si⁺
      0.0,   # 185: SiC2
      0.0,   # 186: SiC2⁺
      0.0,   # 187: SiC2H
      0.0,   # 188: SiC2H⁺
      0.0,   # 189: SiC3
      0.0,   # 190: SiC3⁺
      0.0,   # 191: SiC
      0.0,   # 192: SiC⁺
      0.0,   # 193: SiCH2
      0.0,   # 194: SiCH2⁺
      0.0,   # 195: SiCH3
      0.0,   # 196: SiCH3⁺
      0.0,   # 197: SiH2
      0.0,   # 198: SiH2⁺
      0.0,   # 199: SiH3
      0.0,   # 200: SiH3⁺
      0.0,   # 201: SiN
      0.0,   # 202: SiN⁺
      0.0,   # 203: SiS
      0.0,   # 204: SiS⁺
      0.0,   # 205: C2⁺
      0.0,   # 206: S
      0.0,   # 207: S⁺
      0.0,   # 208: CN⁺
      0.0,   # 209: CO⁺
      0.0,   # 210: N2⁺
      0.0,   # 211: N2
      0.0,   # 212: O2⁺
      0.0,   # 213: C2H⁺
      0.0,   # 214: C2H2⁺
      0.0,   # 215: C2H3
      0.0,   # 216: C2H3⁺
      0.0,   # 217: C5H2⁺
      0.0,   # 218: C6H2⁺
      0.0,   # 219: C7H2⁺
      0.0,   # 220: C3H3⁺
      0.0,   # 221: CO2⁺
      0.0,   # 222: HC3N⁺
      0.0,   # 223: HCN⁺
      0.0,   # 224: C2N2⁺
      0.0,   # 225: C3⁺
      0.0,   # 226: C5⁺
      0.0,   # 227: CH3CH3⁺
      0.0,   # 228: CH3CH3
      0.0,   # 229: PN⁺
      0.0,   # 230: PN
      0.0,   # 231: H2O⁺
      0.0,   # 232: NH2⁺
      0.0,   # 233: O⁺
      0.0,   # 234: OH⁺
      0.0,   # 235: CH3⁺
      0.0,   # 236: CH4⁺
      0.0,   # 237: CH3OH⁺
      0.0,   # 238: N⁺
      0.0,   # 239: SO2⁺
      0.0,   # 240: CS⁺
      0.0,   # 241: Cl⁺
      0.0,   # 242: Cl
      0.0,   # 243: C10⁺
      0.0,   # 244: C10H2⁺
      0.0,   # 245: C3H2
      0.0,   # 246: C3H2⁺
      0.0,   # 247: C3H⁺
      0.0,   # 248: C4⁺
      0.0,   # 249: C4H⁺
      0.0,   # 250: C4P
      0.0,   # 251: C4P⁺
      0.0,   # 252: C5H⁺
      0.0,   # 253: C6⁺
      0.0,   # 254: C6H⁺
      0.0,   # 255: C7⁺
      0.0,   # 256: C7H⁺
      0.0,   # 257: C8⁺
      0.0,   # 258: C8H2⁺
      0.0,   # 259: C8H⁺
      0.0,   # 260: C9⁺
      0.0,   # 261: C9H2⁺
      0.0,   # 262: C9H⁺
      0.0,   # 263: CH3C4H
      0.0,   # 264: CH3C4H⁺
      0.0,   # 265: CH3C6H
      0.0,   # 266: C7H4⁺
      0.0,   # 267: CH3CN⁺
      0.0,   # 268: H2CS
      0.0,   # 269: H2CS⁺
      0.0,   # 270: H2S2
      0.0,   # 271: H2S2⁺
      0.0,   # 272: HC2P
      0.0,   # 273: HC2P⁺
      0.0,   # 274: HC5N⁺
      0.0,   # 275: HC7N⁺
      0.0,   # 276: HC9N⁺
      0.0,   # 277: HCSi
      0.0,   # 278: HCSi⁺
      0.0,   # 279: HCl
      0.0,   # 280: HCl⁺
      0.0,   # 281: HNSi
      0.0,   # 282: HNSi⁺
      0.0,   # 283: HS2
      0.0,   # 284: HS2⁺
      0.0,   # 285: HS⁺
      0.0,   # 286: N2O
      0.0,   # 287: N2O⁺
      0.0,   # 288: NH⁺
      0.0,   # 289: PH2
      0.0,   # 290: PH2⁺
      0.0,   # 291: S2
      0.0,   # 292: S2⁺
      0.0,   # 293: SiC2H2
      0.0,   # 294: SiC2H2⁺
      0.0,   # 295: SiC3H
      0.0,   # 296: SiC3H⁺
      0.0,   # 297: SiC4
      0.0,   # 298: SiC4⁺
      0.0,   # 299: SiH4
      0.0,   # 300: SiH4⁺
      0.0,   # 301: SiH
      0.0,   # 302: SiH⁺
      0.0,   # 303: SiNC
      0.0,   # 304: SiNC⁺
      0.0,   # 305: SiO
      0.0,   # 306: SiO⁺
      0.0,   # 307: H2⁺
      0.0,   # 308: F⁺
      0.0,   # 309: F
      0.0,   # 310: He⁺
      0.0,   # 311: He
      0.0,   # 312: C3N⁺
      0.0,   # 313: HCNO
      0.0,   # 314: HCNO⁺
      0.0,   # 315: HNCO
      0.0,   # 316: HNCO⁺
      0.0,   # 317: HONC
      0.0,   # 318: HONC⁺
      0.0,   # 319: HNO⁺
      0.0,   # 320: HNO
      0.0,   # 321: HCOOH
      0.0,   # 322: HCOOH⁺
      0.0,   # 323: NO2⁺
      0.0,   # 324: C11
      0.0,   # 325: C2H5CN
      0.0,   # 326: CH2CHCNH⁺
      0.0,   # 327: C3P
      0.0,   # 328: CH2CHCCH
      0.0,   # 329: CH2CHCHCH2
      0.0,   # 330: CH2CHCN
      0.0,   # 331: CH2NH
      0.0,   # 332: CH2PH
      0.0,   # 333: CH3C3N
      0.0,   # 334: CH3C5N
      0.0,   # 335: CH3C7N
      0.0,   # 336: CNO
      0.0,   # 337: H2CN
      0.0,   # 338: H2O2
      0.0,   # 339: HCS
      0.0,   # 340: HCS⁺
      0.0,   # 341: HF
      0.0,   # 342: HNC3
      0.0,   # 343: HOCN
      0.0,   # 344: NH2CN
      0.0,   # 345: O2H
      0.0,   # 346: OCN
      0.0,   # 347: SiO2
      0.0,   # 348: C10H3⁺
      0.0]   # 349: C11⁺


### Parameters ###
T = 10
omega = 0.5
params = Dict(
         :T => T,
         :Av => 2,
         :n_H => 611,
         :cr_ion_rate => 6e-18,
         :omega => omega,
         :k1 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 1 with type AD
         :k2 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 2 with type AD
         :k3 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 3 with type AD
         :k4 => 4.7e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 4 with type AD
         :k5 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 5 with type AD
         :k6 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 6 with type AD
         :k7 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 7 with type AD
         :k8 => 5e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 8 with type AD
         :k9 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 9 with type AD
         :k10 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 10 with type AD
         :k11 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 11 with type AD
         :k12 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 12 with type AD
         :k13 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 13 with type AD
         :k14 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 14 with type AD
         :k15 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 15 with type AD
         :k16 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 16 with type AD
         :k17 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 17 with type AD
         :k18 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 18 with type AD
         :k19 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 19 with type AD
         :k20 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 20 with type AD
         :k21 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 21 with type AD
         :k22 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 22 with type AD
         :k23 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 23 with type AD
         :k24 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 24 with type AD
         :k25 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 25 with type AD
         :k26 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 26 with type AD
         :k27 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 27 with type AD
         :k28 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 28 with type AD
         :k29 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 29 with type AD
         :k30 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 30 with type AD
         :k31 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 31 with type AD
         :k32 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 32 with type AD
         :k33 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 33 with type AD
         :k34 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 34 with type AD
         :k35 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 35 with type AD
         :k36 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 36 with type AD
         :k37 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 37 with type AD
         :k38 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 38 with type AD
         :k39 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 39 with type AD
         :k40 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 40 with type AD
         :k41 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 41 with type AD
         :k42 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 42 with type AD
         :k43 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 43 with type AD
         :k44 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 44 with type AD
         :k45 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 45 with type AD
         :k46 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 46 with type AD
         :k47 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 47 with type AD
         :k48 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 48 with type AD
         :k49 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 49 with type AD
         :k50 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 50 with type AD
         :k51 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 51 with type AD
         :k52 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 52 with type AD
         :k53 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 53 with type AD
         :k54 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 54 with type AD
         :k55 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 55 with type AD
         :k56 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 56 with type AD
         :k57 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 57 with type AD
         :k58 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 58 with type AD
         :k59 => 1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 59 with type AD
         :k60 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 60 with type AD
         :k61 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 61 with type AD
         :k62 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 62 with type AD
         :k63 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 63 with type AD
         :k64 => 1.09e-11 * (T/300)^(-2.19) * exp(-(165.1/T)),  # Reaction rate number 64 with type AD
         :k65 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 65 with type AD
         :k66 => 3e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 66 with type AD
         :k67 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 67 with type AD
         :k68 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 68 with type AD
         :k69 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 69 with type AD
         :k70 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 70 with type AD
         :k71 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 71 with type AD
         :k72 => 1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 72 with type AD
         :k73 => 1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 73 with type AD
         :k74 => 2e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 74 with type AD
         :k75 => 4.82e-09 * (T/300)^(0.02) * exp(-(4.3/T)),  # Reaction rate number 75 with type AD
         :k76 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 76 with type AD
         :k77 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 77 with type AD
         :k78 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 78 with type AD
         :k79 => 1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 79 with type AD
         :k80 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 80 with type AD
         :k81 => 1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 81 with type AD
         :k82 => 1e-13 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 82 with type AD
         :k83 => 7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 83 with type AD
         :k84 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 84 with type AD
         :k85 => 1.8e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 85 with type AD
         :k86 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 86 with type AD
         :k87 => 7.7e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 87 with type AD
         :k88 => 1.6e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 88 with type AD
         :k89 => 7.7e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 89 with type AD
         :k90 => 1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 90 with type AD
         :k91 => 5.4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 91 with type AD
         :k92 => 6.2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 92 with type AD
         :k93 => 8.3e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 93 with type AD
         :k94 => 6.2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 94 with type AD
         :k95 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 95 with type AD
         :k96 => 5.8e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 96 with type AD
         :k97 => 6.1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 97 with type AD
         :k98 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 98 with type AD
         :k99 => 2.83e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 99 with type AD
         :k100 => 7.4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 100 with type AD
         :k101 => 2.41e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 101 with type AD
         :k102 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 102 with type AD
         :k103 => 1.22e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 103 with type AD
         :k104 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 104 with type AD
         :k105 => 1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 105 with type AD
         :k106 => 6.3e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 106 with type AD
         :k107 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 107 with type AD
         :k108 => 1.4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 108 with type AD
         :k109 => 1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 109 with type AD
         :k110 => 1.15e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 110 with type AD
         :k111 => 3e-12 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 111 with type AD
         :k112 => 1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 112 with type AD
         :k113 => 5e-12 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 113 with type AD
         :k114 => 1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 114 with type AD
         :k115 => 3e-12 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 115 with type AD
         :k116 => 1.35e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 116 with type AD
         :k117 => 5e-12 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 117 with type AD
         :k118 => 1.1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 118 with type AD
         :k119 => 5e-12 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 119 with type AD
         :k120 => 1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 120 with type AD
         :k121 => 5e-12 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 121 with type AD
         :k122 => 2.2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 122 with type AD
         :k123 => 1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 123 with type AD
         :k124 => 6.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 124 with type AD
         :k125 => 3.1e-10 * (T/300)^(-0.83) * exp(-(0.0/T)),  # Reaction rate number 125 with type AD
         :k126 => 1.9e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 126 with type AD
         :k127 => 3e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 127 with type AD
         :k128 => 2.9e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 128 with type AD
         :k129 => 3.1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 129 with type AD
         :k130 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 130 with type AD
         :k131 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 131 with type AD
         :k132 => 1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 132 with type AD
         :k133 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 133 with type CD
         :k134 => 6e-09 * (T/300)^(0.0) * exp(-(40200.0/T)),  # Reaction rate number 134 with type CD
         :k135 => 1e-08 * (T/300)^(0.0) * exp(-(84100.0/T)),  # Reaction rate number 135 with type CD
         :k136 => 5.8e-09 * (T/300)^(0.0) * exp(-(52900.0/T)),  # Reaction rate number 136 with type CD
         :k137 => 3.8e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 137 with type CD
         :k138 => 6e-09 * (T/300)^(0.0) * exp(-(52300.0/T)),  # Reaction rate number 138 with type CD
         :k139 => 6e-09 * (T/300)^(0.0) * exp(-(50900.0/T)),  # Reaction rate number 139 with type CD
         :k140 => 3.22e-09 * (T/300)^(0.35) * exp(-(102000.0/T)),  # Reaction rate number 140 with type CD
         :k141 => 6e-09 * (T/300)^(0.0) * exp(-(40200.0/T)),  # Reaction rate number 141 with type CD
         :k142 => 4.67e-07 * (T/300)^(-1.0) * exp(-(55000.0/T)),  # Reaction rate number 142 with type CD
         :k143 => 5.8e-09 * (T/300)^(0.0) * exp(-(52900.0/T)),  # Reaction rate number 143 with type CD
         :k144 => 6e-09 * (T/300)^(0.0) * exp(-(52300.0/T)),  # Reaction rate number 144 with type CD
         :k145 => 6e-09 * (T/300)^(0.0) * exp(-(50900.0/T)),  # Reaction rate number 145 with type CD
         :k146 => 7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 146 with type CD
         :k147 => 5e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 147 with type CE
         :k148 => 1.7e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 148 with type CE
         :k149 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 149 with type CE
         :k150 => 1.5e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 150 with type CE
         :k151 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 151 with type CE
         :k152 => 5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 152 with type CE
         :k153 => 3.33e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 153 with type CE
         :k154 => 5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 154 with type CE
         :k155 => 3.1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 155 with type CE
         :k156 => 5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 156 with type CE
         :k157 => 1.61e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 157 with type CE
         :k158 => 5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 158 with type CE
         :k159 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 159 with type CE
         :k160 => 5.2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 160 with type CE
         :k161 => 2.1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 161 with type CE
         :k162 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 162 with type CE
         :k163 => 2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 163 with type CE
         :k164 => 2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 164 with type CE
         :k165 => 5.7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 165 with type CE
         :k166 => 2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 166 with type CE
         :k167 => 1.5e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 167 with type CE
         :k168 => 4.15e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 168 with type CE
         :k169 => 2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 169 with type CE
         :k170 => 3.8e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 170 with type CE
         :k171 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 171 with type CE
         :k172 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 172 with type CE
         :k173 => 2.6e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 173 with type CE
         :k174 => 7.8e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 174 with type CE
         :k175 => 6e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 175 with type CE
         :k176 => 2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 176 with type CE
         :k177 => 1.31e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 177 with type CE
         :k178 => 4.8e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 178 with type CE
         :k179 => 3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 179 with type CE
         :k180 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 180 with type CE
         :k181 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 181 with type CE
         :k182 => 1.1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 182 with type CE
         :k183 => 2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 183 with type CE
         :k184 => 1.7e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 184 with type CE
         :k185 => 6.72e-10 * (T/300)^(0.0) * exp(-(-0.5/T)),  # Reaction rate number 185 with type CE
         :k186 => 7.05e-10 * (T/300)^(-0.03) * exp(-(-16.7/T)),  # Reaction rate number 186 with type CE
         :k187 => 7.6e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 187 with type CE
         :k188 => 1.1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 188 with type CE
         :k189 => 4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 189 with type CE
         :k190 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 190 with type CE
         :k191 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 191 with type CE
         :k192 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 192 with type CE
         :k193 => 2.6e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 193 with type CE
         :k194 => 2.1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 194 with type CE
         :k195 => 2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 195 with type CE
         :k196 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 196 with type CE
         :k197 => 2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 197 with type CE
         :k198 => 2.5e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 198 with type CE
         :k199 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 199 with type CE
         :k200 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 200 with type CE
         :k201 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 201 with type CE
         :k202 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 202 with type CE
         :k203 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 203 with type CE
         :k204 => 2.3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 204 with type CE
         :k205 => 3.8e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 205 with type CE
         :k206 => 3.4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 206 with type CE
         :k207 => 5.8e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 207 with type CE
         :k208 => 8.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 208 with type CE
         :k209 => 8.4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 209 with type CE
         :k210 => 8.4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 210 with type CE
         :k211 => 4.1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 211 with type CE
         :k212 => 1.2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 212 with type CE
         :k213 => 1.2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 213 with type CE
         :k214 => 3.3e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 214 with type CE
         :k215 => 4.14e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 215 with type CE
         :k216 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 216 with type CE
         :k217 => 5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 217 with type CE
         :k218 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 218 with type CE
         :k219 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 219 with type CE
         :k220 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 220 with type CE
         :k221 => 2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 221 with type CE
         :k222 => 8.6e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 222 with type CE
         :k223 => 2.2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 223 with type CE
         :k224 => 1.26e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 224 with type CE
         :k225 => 5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 225 with type CE
         :k226 => 1.2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 226 with type CE
         :k227 => 2.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 227 with type CE
         :k228 => 7.3e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 228 with type CE
         :k229 => 1.28e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 229 with type CE
         :k230 => 1.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 230 with type CE
         :k231 => 5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 231 with type CE
         :k232 => 8e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 232 with type CE
         :k233 => 6.57e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 233 with type CE
         :k234 => 1.3e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 234 with type CE
         :k235 => 3.96e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 235 with type CE
         :k236 => 3.06e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 236 with type CE
         :k237 => 1.15e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 237 with type CE
         :k238 => 1.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 238 with type CE
         :k239 => 5.36e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 239 with type CE
         :k240 => 6.8e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 240 with type CE
         :k241 => 8e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 241 with type CE
         :k242 => 3.9e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 242 with type CE
         :k243 => 7.9e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 243 with type CE
         :k244 => 1.1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 244 with type CE
         :k245 => 1.1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 245 with type CE
         :k246 => 1.1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 246 with type CE
         :k247 => 1.1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 247 with type CE
         :k248 => 5.2e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 248 with type CE
         :k249 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 249 with type CE
         :k250 => 2.6e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 250 with type CE
         :k251 => 4.6e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 251 with type CE
         :k252 => 3.6e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 252 with type CE
         :k253 => 4.59e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 253 with type CE
         :k254 => 7.6e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 254 with type CE
         :k255 => 3.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 255 with type CE
         :k256 => 4.7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 256 with type CE
         :k257 => 2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 257 with type CE
         :k258 => 4.2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 258 with type CE
         :k259 => 4.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 259 with type CE
         :k260 => 8.8e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 260 with type CE
         :k261 => 4.3e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 261 with type CE
         :k262 => 4.3e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 262 with type CE
         :k263 => 4.7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 263 with type CE
         :k264 => 8.7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 264 with type CE
         :k265 => 4.9e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 265 with type CE
         :k266 => 9.7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 266 with type CE
         :k267 => 4.3e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 267 with type CE
         :k268 => 4.8e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 268 with type CE
         :k269 => 3e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 269 with type CE
         :k270 => 2.4e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 270 with type CE
         :k271 => 2.4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 271 with type CE
         :k272 => 4.4e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 272 with type CE
         :k273 => 3.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 273 with type CE
         :k274 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 274 with type CE
         :k275 => 3.4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 275 with type CE
         :k276 => 1.98e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 276 with type CE
         :k277 => 1.13e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 277 with type CE
         :k278 => 1.38e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 278 with type CE
         :k279 => 1.8e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 279 with type CE
         :k280 => 1.62e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 280 with type CE
         :k281 => 9.45e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 281 with type CE
         :k282 => 1.65e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 282 with type CE
         :k283 => 3.9e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 283 with type CE
         :k284 => 4.2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 284 with type CE
         :k285 => 7.93e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 285 with type CE
         :k286 => 5.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 286 with type CE
         :k287 => 3.2e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 287 with type CE
         :k288 => 6.4e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 288 with type CE
         :k289 => 3.2e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 289 with type CE
         :k290 => 3.1e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 290 with type CE
         :k291 => 3.4e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 291 with type CE
         :k292 => 3.6e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 292 with type CE
         :k293 => 6.3e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 293 with type CE
         :k294 => 3.5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 294 with type CE
         :k295 => 3.5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 295 with type CE
         :k296 => 3.1e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 296 with type CE
         :k297 => 3.5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 297 with type CE
         :k298 => 3e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 298 with type CE
         :k299 => 6.3e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 299 with type CE
         :k300 => 5.2e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 300 with type CE
         :k301 => 1.79e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 301 with type CE
         :k302 => 3.7e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 302 with type CE
         :k303 => 5.7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 303 with type CE
         :k304 => 2.58e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 304 with type CE
         :k305 => 1.1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 305 with type CE
         :k306 => 1e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 306 with type CE
         :k307 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 307 with type CE
         :k308 => 1.35e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 308 with type CE
         :k309 => 2.44e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 309 with type CE
         :k310 => 7.4e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 310 with type CE
         :k311 => 3.3e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 311 with type CE
         :k312 => 1.2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 312 with type CE
         :k313 => 1.1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 313 with type CE
         :k314 => 9.6e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 314 with type CE
         :k315 => 1.5e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 315 with type CE
         :k316 => 7.4e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 316 with type CE
         :k317 => 1.7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 317 with type CE
         :k318 => 6.2e-11 * (T/300)^(0.79) * exp(-(6920.0/T)),  # Reaction rate number 318 with type CE
         :k319 => 9.3e-11 * (T/300)^(0.73) * exp(-(232.0/T)),  # Reaction rate number 319 with type CE
         :k320 => 2.76e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 320 with type CE
         :k321 => 2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 321 with type CE
         :k322 => 3.98e-08 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 322 with type CE
         :k323 => 3.1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 323 with type CE
         :k324 => 5.4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 324 with type CE
         :k325 => 2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 325 with type CE
         :k326 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 326 with type CE
         :k327 => 3.06e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 327 with type CE
         :k328 => 3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 328 with type CE
         :k329 => 1.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 329 with type CE
         :k330 => 4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 330 with type CE
         :k331 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 331 with type CE
         :k332 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 332 with type CE
         :k333 => 4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 333 with type CE
         :k334 => 2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 334 with type CE
         :k335 => 2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 335 with type CE
         :k336 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 336 with type CE
         :k337 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 337 with type CE
         :k338 => 4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 338 with type CE
         :k339 => 2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 339 with type CE
         :k340 => 2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 340 with type CE
         :k341 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 341 with type CE
         :k342 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 342 with type CE
         :k343 => 4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 343 with type CE
         :k344 => 2.5e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 344 with type CE
         :k345 => 3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 345 with type CE
         :k346 => 4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 346 with type CE
         :k347 => 2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 347 with type CE
         :k348 => 2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 348 with type CE
         :k349 => 4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 349 with type CE
         :k350 => 2.5e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 350 with type CE
         :k351 => 3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 351 with type CE
         :k352 => 4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 352 with type CE
         :k353 => 2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 353 with type CE
         :k354 => 2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 354 with type CE
         :k355 => 4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 355 with type CE
         :k356 => 2.5e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 356 with type CE
         :k357 => 3.5e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 357 with type CE
         :k358 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 358 with type CE
         :k359 => 1.4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 359 with type CE
         :k360 => 2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 360 with type CE
         :k361 => 6.3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 361 with type CE
         :k362 => 3.4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 362 with type CE
         :k363 => 2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 363 with type CE
         :k364 => 2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 364 with type CE
         :k365 => 2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 365 with type CE
         :k366 => 3.6e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 366 with type CE
         :k367 => 8.4e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 367 with type CE
         :k368 => 4e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 368 with type CE
         :k369 => 2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 369 with type CE
         :k370 => 5.9e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 370 with type CE
         :k371 => 1.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 371 with type CE
         :k372 => 1.9e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 372 with type CE
         :k373 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 373 with type CE
         :k374 => 4.9e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 374 with type CE
         :k375 => 7.4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 375 with type CE
         :k376 => 2.96e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 376 with type CE
         :k377 => 4.7e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 377 with type CE
         :k378 => 6.9e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 378 with type CE
         :k379 => 3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 379 with type CE
         :k380 => 5.28e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 380 with type CE
         :k381 => 1.5e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 381 with type CE
         :k382 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 382 with type CE
         :k383 => 4e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 383 with type CE
         :k384 => 2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 384 with type CE
         :k385 => 4e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 385 with type CE
         :k386 => 1e-08 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 386 with type CE
         :k387 => 1e-08 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 387 with type CE
         :k388 => 1.05e-08 * (T/300)^(-0.13) * exp(-(0.0/T)),  # Reaction rate number 388 with type CE
         :k389 => 9.4e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 389 with type CE
         :k390 => 3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 390 with type CE
         :k391 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 391 with type CE
         :k392 => 1.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 392 with type CE
         :k393 => 3.3e-09 * (T/300)^(1.0) * exp(-(0.0/T)),  # Reaction rate number 393 with type CE
         :k394 => 1.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 394 with type CE
         :k395 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 395 with type CE
         :k396 => 3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 396 with type CE
         :k397 => 1.6e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 397 with type CE
         :k398 => 1.1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 398 with type CE
         :k399 => 1.85e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 399 with type CE
         :k400 => 2.9e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 400 with type CE
         :k401 => 3.7e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 401 with type CE
         :k402 => 2.1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 402 with type CE
         :k403 => 2.9e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 403 with type CE
         :k404 => 4.7e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 404 with type CE
         :k405 => 2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 405 with type CE
         :k406 => 6.86e-10 * (T/300)^(0.26) * exp(-(224.3/T)),  # Reaction rate number 406 with type CE
         :k407 => 2.1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 407 with type CE
         :k408 => 2.1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 408 with type CE
         :k409 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 409 with type CE
         :k410 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 410 with type CE
         :k411 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 411 with type CE
         :k412 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 412 with type CE
         :k413 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 413 with type CE
         :k414 => 3e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 414 with type CE
         :k415 => 1.3e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 415 with type CE
         :k416 => 5.78e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 416 with type CE
         :k417 => 3.2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 417 with type CE
         :k418 => 9.9e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 418 with type CE
         :k419 => 3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 419 with type CE
         :k420 => 1.5e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 420 with type CE
         :k421 => 1.5e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 421 with type CE
         :k422 => 3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 422 with type CE
         :k423 => 1.5e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 423 with type CE
         :k424 => 1.5e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 424 with type CE
         :k425 => 3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 425 with type CE
         :k426 => 1.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 426 with type CE
         :k427 => 1.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 427 with type CE
         :k428 => 1.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 428 with type CE
         :k429 => 1.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 429 with type CE
         :k430 => 1.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 430 with type CE
         :k431 => 1.7e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 431 with type CE
         :k432 => 3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 432 with type CE
         :k433 => 3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 433 with type CE
         :k434 => 3.3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 434 with type CE
         :k435 => 1.5e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 435 with type CE
         :k436 => 1.1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 436 with type CE
         :k437 => 4.82e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 437 with type CE
         :k438 => 2.21e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 438 with type CE
         :k439 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 439 with type CE
         :k440 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 440 with type CE
         :k441 => 2.94e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 441 with type CE
         :k442 => 1.4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 442 with type CE
         :k443 => 7.1e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 443 with type CE
         :k444 => 1.2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 444 with type CE
         :k445 => 6.44e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 445 with type CE
         :k446 => 1.4e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 446 with type CE
         :k447 => 3.9e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 447 with type CE
         :k448 => 2.7e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 448 with type CE
         :k449 => 2.7e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 449 with type CE
         :k450 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 450 with type CE
         :k451 => 2.1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 451 with type CE
         :k452 => 5.7e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 452 with type CE
         :k453 => 7.6e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 453 with type CE
         :k454 => 1.1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 454 with type CE
         :k455 => 8e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 455 with type CE
         :k456 => 7.6e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 456 with type CE
         :k457 => 6.24e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 457 with type CE
         :k458 => 7.2e-15 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 458 with type CE
         :k459 => 1.9e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 459 with type CE
         :k460 => 5.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 460 with type CE
         :k461 => 2.07e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 461 with type CE
         :k462 => 4.7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 462 with type CE
         :k463 => 1.9e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 463 with type CE
         :k464 => 1.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 464 with type CE
         :k465 => 4.4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 465 with type CE
         :k466 => 6.4e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 466 with type CE
         :k467 => 1.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 467 with type CE
         :k468 => 1.41e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 468 with type CE
         :k469 => 9.72e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 469 with type CE
         :k470 => 2.8e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 470 with type CE
         :k471 => 2.2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 471 with type CE
         :k472 => 2.7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 472 with type CE
         :k473 => 6.2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 473 with type CE
         :k474 => 4.6e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 474 with type CE
         :k475 => 2.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 475 with type CE
         :k476 => 4.3e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 476 with type CE
         :k477 => 3e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 477 with type CE
         :k478 => 1.72e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 478 with type CE
         :k479 => 2.04e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 479 with type CE
         :k480 => 1.8e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 480 with type CE
         :k481 => 2.3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 481 with type CE
         :k482 => 1.89e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 482 with type CE
         :k483 => 1.8e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 483 with type CE
         :k484 => 1.4e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 484 with type CE
         :k485 => 8.45e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 485 with type CE
         :k486 => 1.4e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 486 with type CE
         :k487 => 6.4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 487 with type CE
         :k488 => 7.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 488 with type CE
         :k489 => 6.4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 489 with type CE
         :k490 => 3.7e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 490 with type CE
         :k491 => 1.2e-15 * (T/300)^(0.25) * exp(-(0.0/T)),  # Reaction rate number 491 with type CE
         :k492 => 5.66e-10 * (T/300)^(0.36) * exp(-(-8.6/T)),  # Reaction rate number 492 with type CE
         :k493 => 8.9e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 493 with type CE
         :k494 => 8.1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 494 with type CE
         :k495 => 3.2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 495 with type CE
         :k496 => 5.7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 496 with type CE
         :k497 => 3.4e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 497 with type CE
         :k498 => 3.9e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 498 with type CE
         :k499 => 1.9e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 499 with type CE
         :k500 => 3.6e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 500 with type CE
         :k501 => 7e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 501 with type CE
         :k502 => 3.6e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 502 with type CE
         :k503 => 3.6e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 503 with type CE
         :k504 => 6.6e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 504 with type CE
         :k505 => 1.6e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 505 with type CE
         :k506 => 1.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 506 with type CE
         :k507 => 1.5e-08 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 507 with type CE
         :k508 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 508 with type CE
         :k509 => 2.54e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 509 with type CE
         :k510 => 2.4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 510 with type CE
         :k511 => 5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 511 with type CE
         :k512 => 8.2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 512 with type CE
         :k513 => 6.3e-15 * (T/300)^(0.75) * exp(-(0.0/T)),  # Reaction rate number 513 with type CE
         :k514 => 5.1e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 514 with type CE
         :k515 => 5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 515 with type CE
         :k516 => 1.21e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 516 with type CE
         :k517 => 9.69e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 517 with type CE
         :k518 => 6.05e-11 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 518 with type CE
         :k519 => 3.08e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 519 with type CE
         :k520 => 8.39e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 520 with type CE
         :k521 => 5.68e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 521 with type CE
         :k522 => 1.03e-08 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 522 with type CE
         :k523 => 6.4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 523 with type CE
         :k524 => 2.64e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 524 with type CE
         :k525 => 3.3e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 525 with type CE
         :k526 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 526 with type CE
         :k527 => 4.3e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 527 with type CE
         :k528 => 3.3e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 528 with type CE
         :k529 => 3e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 529 with type CE
         :k530 => 2.7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 530 with type CE
         :k531 => 2.9e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 531 with type CE
         :k532 => 2.8e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 532 with type CE
         :k533 => 2.9e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 533 with type CE
         :k534 => 2.6e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 534 with type CE
         :k535 => 7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 535 with type CE
         :k536 => 8.1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 536 with type CE
         :k537 => 1.2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 537 with type CE
         :k538 => 2.8e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 538 with type CE
         :k539 => 1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 539 with type CE
         :k540 => 2.9e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 540 with type CE
         :k541 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 541 with type CE
         :k542 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 542 with type CE
         :k543 => 9.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 543 with type CE
         :k544 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 544 with type CE
         :k545 => 1.24e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 545 with type CE
         :k546 => 2.8e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 546 with type CE
         :k547 => 1.1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 547 with type CE
         :k548 => 7.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 548 with type CE
         :k549 => 8.25e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 549 with type CE
         :k550 => 1.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 550 with type CE
         :k551 => 1.88e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 551 with type CE
         :k552 => 2.8e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 552 with type CE
         :k553 => 1.06e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 553 with type CE
         :k554 => 3.7e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 554 with type CE
         :k555 => 4.5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 555 with type CE
         :k556 => 1.2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 556 with type CE
         :k557 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 557 with type CE
         :k558 => 1.97e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 558 with type CE
         :k559 => 3.7e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 559 with type CE
         :k560 => 4.51e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 560 with type CE
         :k561 => 3.11e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 561 with type CE
         :k562 => 1.02e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 562 with type CE
         :k563 => 3.7e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 563 with type CE
         :k564 => 7.7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 564 with type CE
         :k565 => 4.3e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 565 with type CE
         :k566 => 3.77e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 566 with type CE
         :k567 => 1.5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 567 with type CE
         :k568 => 3.7e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 568 with type CE
         :k569 => 4.4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 569 with type CE
         :k570 => 5e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 570 with type CE
         :k571 => 2.6e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 571 with type CE
         :k572 => 1.1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 572 with type CE
         :k573 => 1e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 573 with type CE
         :k574 => 9.9e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 574 with type CE
         :k575 => 1.05e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 575 with type CE
         :k576 => 1.8e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 576 with type CE
         :k577 => 7.12e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 577 with type CE
         :k578 => 4.51e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 578 with type CE
         :k579 => 6.9e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 579 with type CE
         :k580 => 7.2e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 580 with type CE
         :k581 => 4.3e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 581 with type CE
         :k582 => 6.9e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 582 with type CE
         :k583 => 7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 583 with type CE
         :k584 => 4.4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 584 with type CE
         :k585 => 4.6e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 585 with type CE
         :k586 => 9.1e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 586 with type CE
         :k587 => 4.5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 587 with type CE
         :k588 => 4.9e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 588 with type CE
         :k589 => 8.9e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 589 with type CE
         :k590 => 8.7e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 590 with type CE
         :k591 => 5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 591 with type CE
         :k592 => 2.3e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 592 with type CE
         :k593 => 4.2e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 593 with type CE
         :k594 => 3.3e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 594 with type CE
         :k595 => 7.2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 595 with type CE
         :k596 => 3.2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 596 with type CE
         :k597 => 1.9e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 597 with type CE
         :k598 => 2.14e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 598 with type CE
         :k599 => 1.8e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 599 with type CE
         :k600 => 3.2e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 600 with type CE
         :k601 => 6.24e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 601 with type CE
         :k602 => 2.02e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 602 with type CE
         :k603 => 1.9e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 603 with type CE
         :k604 => 4.25e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 604 with type CE
         :k605 => 2.21e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 605 with type CE
         :k606 => 3.4e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 606 with type CE
         :k607 => 1.7e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 607 with type CE
         :k608 => 1.68e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 608 with type CE
         :k609 => 5.25e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 609 with type CE
         :k610 => 1.9e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 610 with type CE
         :k611 => 2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 611 with type CE
         :k612 => 2.3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 612 with type CE
         :k613 => 3.08e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 613 with type CE
         :k614 => 1.44e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 614 with type CE
         :k615 => 1.3e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 615 with type CE
         :k616 => 6.5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 616 with type CE
         :k617 => 3.2e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 617 with type CE
         :k618 => 6.5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 618 with type CE
         :k619 => 3.6e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 619 with type CE
         :k620 => 9.2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 620 with type CE
         :k621 => 1.3e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 621 with type CE
         :k622 => 3.8e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 622 with type CE
         :k623 => 1.2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 623 with type CE
         :k624 => 7.8e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 624 with type CE
         :k625 => 3.7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 625 with type CE
         :k626 => 7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 626 with type CE
         :k627 => 4.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 627 with type CE
         :k628 => 4.6e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 628 with type CE
         :k629 => 3.7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 629 with type CE
         :k630 => 5.1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 630 with type CE
         :k631 => 7.2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 631 with type CE
         :k632 => 2.7e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 632 with type CE
         :k633 => 2.3e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 633 with type CE
         :k634 => 1e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 634 with type CE
         :k635 => 2.6e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 635 with type CE
         :k636 => 2.5e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 636 with type CE
         :k637 => 2.6e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 637 with type CE
         :k638 => 2.2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 638 with type CE
         :k639 => 1e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 639 with type CE
         :k640 => 2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 640 with type CE
         :k641 => 7.7e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 641 with type CE
         :k642 => 7.1e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 642 with type CE
         :k643 => 2.6e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 643 with type CE
         :k644 => 2.3e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 644 with type CE
         :k645 => 2.7e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 645 with type CE
         :k646 => 4.8e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 646 with type CE
         :k647 => 3.9e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 647 with type CE
         :k648 => 7e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 648 with type CE
         :k649 => 4.6e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 649 with type CE
         :k650 => 2.94e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 650 with type CE
         :k651 => 4.75e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 651 with type CE
         :k652 => 8.9e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 652 with type CE
         :k653 => 4.9e-12 * (T/300)^(0.5) * exp(-(4580.0/T)),  # Reaction rate number 653 with type CE
         :k654 => 2.9e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 654 with type CE
         :k655 => 2.1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 655 with type CE
         :k656 => 3.2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 656 with type CE
         :k657 => 1.36e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 657 with type CE
         :k658 => 4.3e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 658 with type CE
         :k659 => 6.3e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 659 with type CE
         :k660 => 1e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 660 with type CE
         :k661 => 1.2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 661 with type CE
         :k662 => 1.9e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 662 with type CE
         :k663 => 6.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 663 with type CE
         :k664 => 3.6e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 664 with type CE
         :k665 => 2.04e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 665 with type CE
         :k666 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 666 with type CE
         :k667 => 7.3e-10 * (T/300)^(0.0) * exp(-(890.0/T)),  # Reaction rate number 667 with type CE
         :k668 => 1.11e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 668 with type CE
         :k669 => 1.22e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 669 with type CE
         :k670 => 7.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 670 with type CE
         :k671 => 4.5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 671 with type CE
         :k672 => 5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 672 with type CE
         :k673 => 1.1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 673 with type CE
         :k674 => 1.4e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 674 with type CE
         :k675 => 1.26e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 675 with type CE
         :k676 => 6.6e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 676 with type CE
         :k677 => 5.4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 677 with type CE
         :k678 => 5.3e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 678 with type CE
         :k679 => 4.6e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 679 with type CE
         :k680 => 2.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 680 with type CE
         :k681 => 6.5e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 681 with type CE
         :k682 => 1.4e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 682 with type CE
         :k683 => 9.62e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 683 with type CE
         :k684 => 1e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 684 with type CE
         :k685 => 4.8e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 685 with type CE
         :k686 => 4.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 686 with type CE
         :k687 => 4.8e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 687 with type CE
         :k688 => 7.44e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 688 with type CE
         :k689 => 1.59e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 689 with type CE
         :k690 => 1.23e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 690 with type CE
         :k691 => 2.8e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 691 with type CE
         :k692 => 1.2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 692 with type CE
         :k693 => 3.59e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 693 with type CE
         :k694 => 5.9e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 694 with type CE
         :k695 => 4.3e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 695 with type CE
         :k696 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 696 with type CE
         :k697 => 6.5e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 697 with type CE
         :k698 => 6.4e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 698 with type CE
         :k699 => 3.1e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 699 with type CE
         :k700 => 6.3e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 700 with type CE
         :k701 => 4.34e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 701 with type CE
         :k702 => 1.44e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 702 with type CE
         :k703 => 1.8e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 703 with type CE
         :k704 => 7.2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 704 with type CE
         :k705 => 3.7e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 705 with type CE
         :k706 => 3.2e-09 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 706 with type CE
         :k707 => 5e-11 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 707 with type CE
         :k708 => 1.1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 708 with type CE
         :k709 => 9.7e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 709 with type CE
         :k710 => 6.84e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 710 with type CE
         :k711 => 1.82e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 711 with type CE
         :k712 => 1.26e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 712 with type CE
         :k713 => 1.8e-10 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 713 with type CE
         :k714 => 1.6e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 714 with type CE
         :k715 => 1.9e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 715 with type CE
         :k716 => 1.5e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 716 with type CE
         :k717 => 2e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 717 with type CE
         :k718 => 1.6e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 718 with type CE
         :k719 => 1.4e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 719 with type CE
         :k720 => 1.6e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 720 with type CE
         :k721 => 1.6e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 721 with type CE
         :k722 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 722 with type CE
         :k723 => 1.6e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 723 with type CE
         :k724 => 4.2e-10 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 724 with type CE
         :k725 => 1e-09 * (T/300)^(0.0) * exp(-(0.0/T)),  # Reaction rate number 725 with type CE
         :k726 => 2.3e-17,  # Reaction rate number 726 with type CP
         :k727 => 3.9e-17,  # Reaction rate number 727 with type CP
         :k728 => 3.9e-17,  # Reaction rate number 728 with type CP
         :k729 => 3.9e-21,  # Reaction rate number 729 with type CP
         :k730 => 2.86e-19,  # Reaction rate number 730 with type CP
         :k731 => 1.2e-17,  # Reaction rate number 731 with type CP
         :k732 => 1.3e-18,  # Reaction rate number 732 with type CP
         :k733 => 5.98e-18,  # Reaction rate number 733 with type CP
         :k734 => 6.5e-18,  # Reaction rate number 734 with type CP
         :k735 => 2.7e-17,  # Reaction rate number 735 with type CP
         :k736 => 3.4e-17,  # Reaction rate number 736 with type CP
         :k737 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 737 with type CR
         :k738 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 738 with type CR
         :k739 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 739 with type CR
         :k740 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 740 with type CR
         :k741 => 1.3e-17 * (T/300)^(0.0) * (875.0/(1-omega)),  # Reaction rate number 741 with type CR
         :k742 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 742 with type CR
         :k743 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 743 with type CR
         :k744 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 744 with type CR
         :k745 => 1.3e-17 * (T/300)^(0.0) * (119.5/(1-omega)),  # Reaction rate number 745 with type CR
         :k746 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 746 with type CR
         :k747 => 1.3e-17 * (T/300)^(0.0) * (654.5/(1-omega)),  # Reaction rate number 747 with type CR
         :k748 => 1.3e-17 * (T/300)^(0.0) * (2577.5/(1-omega)),  # Reaction rate number 748 with type CR
         :k749 => 1.3e-17 * (T/300)^(0.0) * (1500.0/(1-omega)),  # Reaction rate number 749 with type CR
         :k750 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 750 with type CR
         :k751 => 1.3e-17 * (T/300)^(0.0) * (1881.0/(1-omega)),  # Reaction rate number 751 with type CR
         :k752 => 1.3e-17 * (T/300)^(0.0) * (389.0/(1-omega)),  # Reaction rate number 752 with type CR
         :k753 => 1.3e-17 * (T/300)^(0.0) * (1881.0/(1-omega)),  # Reaction rate number 753 with type CR
         :k754 => 1.3e-17 * (T/300)^(0.0) * (389.0/(1-omega)),  # Reaction rate number 754 with type CR
         :k755 => 1.3e-17 * (T/300)^(0.0) * (1122.5/(1-omega)),  # Reaction rate number 755 with type CR
         :k756 => 1.3e-17 * (T/300)^(0.0) * (2388.0/(1-omega)),  # Reaction rate number 756 with type CR
         :k757 => 1.3e-17 * (T/300)^(0.0) * (2153.5/(1-omega)),  # Reaction rate number 757 with type CR
         :k758 => 1.3e-17 * (T/300)^(0.0) * (1368.0/(1-omega)),  # Reaction rate number 758 with type CR
         :k759 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 759 with type CR
         :k760 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 760 with type CR
         :k761 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 761 with type CR
         :k762 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 762 with type CR
         :k763 => 1.3e-17 * (T/300)^(0.0) * (375.0/(1-omega)),  # Reaction rate number 763 with type CR
         :k764 => 1.3e-17 * (T/300)^(0.0) * (375.0/(1-omega)),  # Reaction rate number 764 with type CR
         :k765 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 765 with type CR
         :k766 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 766 with type CR
         :k767 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 767 with type CR
         :k768 => 1.3e-17 * (T/300)^(0.0) * (559.5/(1-omega)),  # Reaction rate number 768 with type CR
         :k769 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 769 with type CR
         :k770 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 770 with type CR
         :k771 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 771 with type CR
         :k772 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 772 with type CR
         :k773 => 1.3e-17 * (T/300)^(0.0) * (875.0/(1-omega)),  # Reaction rate number 773 with type CR
         :k774 => 1.3e-17 * (T/300)^(0.0) * (3304.5/(1-omega)),  # Reaction rate number 774 with type CR
         :k775 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 775 with type CR
         :k776 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 776 with type CR
         :k777 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 777 with type CR
         :k778 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 778 with type CR
         :k779 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 779 with type CR
         :k780 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 780 with type CR
         :k781 => 1.3e-17 * (T/300)^(0.0) * (5000.0/(1-omega)),  # Reaction rate number 781 with type CR
         :k782 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 782 with type CR
         :k783 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 783 with type CR
         :k784 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 784 with type CR
         :k785 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 785 with type CR
         :k786 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 786 with type CR
         :k787 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 787 with type CR
         :k788 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 788 with type CR
         :k789 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 789 with type CR
         :k790 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 790 with type CR
         :k791 => 1.3e-17 * (T/300)^(0.0) * (875.0/(1-omega)),  # Reaction rate number 791 with type CR
         :k792 => 1.3e-17 * (T/300)^(0.0) * (875.0/(1-omega)),  # Reaction rate number 792 with type CR
         :k793 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 793 with type CR
         :k794 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 794 with type CR
         :k795 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 795 with type CR
         :k796 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 796 with type CR
         :k797 => 1.3e-17 * (T/300)^(0.0) * (875.0/(1-omega)),  # Reaction rate number 797 with type CR
         :k798 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 798 with type CR
         :k799 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 799 with type CR
         :k800 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 800 with type CR
         :k801 => 1.3e-17 * (T/300)^(0.0) * (875.0/(1-omega)),  # Reaction rate number 801 with type CR
         :k802 => 1.3e-17 * (T/300)^(0.0) * (1500.0/(1-omega)),  # Reaction rate number 802 with type CR
         :k803 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 803 with type CR
         :k804 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 804 with type CR
         :k805 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 805 with type CR
         :k806 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 806 with type CR
         :k807 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 807 with type CR
         :k808 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 808 with type CR
         :k809 => 1.3e-17 * (T/300)^(0.0) * (875.0/(1-omega)),  # Reaction rate number 809 with type CR
         :k810 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 810 with type CR
         :k811 => 1.3e-17 * (T/300)^(0.0) * (875.0/(1-omega)),  # Reaction rate number 811 with type CR
         :k812 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 812 with type CR
         :k813 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 813 with type CR
         :k814 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 814 with type CR
         :k815 => 1.3e-17 * (T/300)^(0.0) * (875.0/(1-omega)),  # Reaction rate number 815 with type CR
         :k816 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 816 with type CR
         :k817 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 817 with type CR
         :k818 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 818 with type CR
         :k819 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 819 with type CR
         :k820 => 1.3e-17 * (T/300)^(0.0) * (875.0/(1-omega)),  # Reaction rate number 820 with type CR
         :k821 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 821 with type CR
         :k822 => 1.3e-17 * (T/300)^(0.0) * (875.0/(1-omega)),  # Reaction rate number 822 with type CR
         :k823 => 1.3e-17 * (T/300)^(0.0) * (255.0/(1-omega)),  # Reaction rate number 823 with type CR
         :k824 => 1.3e-17 * (T/300)^(0.0) * (375.0/(1-omega)),  # Reaction rate number 824 with type CR
         :k825 => 1.3e-17 * (T/300)^(0.0) * (375.0/(1-omega)),  # Reaction rate number 825 with type CR
         :k826 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 826 with type CR
         :k827 => 1.3e-17 * (T/300)^(0.0) * (88.0/(1-omega)),  # Reaction rate number 827 with type CR
         :k828 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 828 with type CR
         :k829 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 829 with type CR
         :k830 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 830 with type CR
         :k831 => 1.3e-17 * (T/300)^(0.0) * (2652.5/(1-omega)),  # Reaction rate number 831 with type CR
         :k832 => 1.3e-17 * (T/300)^(0.0) * (1642.0/(1-omega)),  # Reaction rate number 832 with type CR
         :k833 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 833 with type CR
         :k834 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 834 with type CR
         :k835 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 835 with type CR
         :k836 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 836 with type CR
         :k837 => 1.3e-17 * (T/300)^(0.0) * (1000.0/(1-omega)),  # Reaction rate number 837 with type CR
         :k838 => 1.3e-17 * (T/300)^(0.0) * (375.0/(1-omega)),  # Reaction rate number 838 with type CR
         :k839 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 839 with type CR
         :k840 => 1.3e-17 * (T/300)^(0.0) * (609.0/(1-omega)),  # Reaction rate number 840 with type CR
         :k841 => 1.3e-17 * (T/300)^(0.0) * (456.5/(1-omega)),  # Reaction rate number 841 with type CR
         :k842 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 842 with type CR
         :k843 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 843 with type CR
         :k844 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 844 with type CR
         :k845 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 845 with type CR
         :k846 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 846 with type CR
         :k847 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 847 with type CR
         :k848 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 848 with type CR
         :k849 => 1.3e-17 * (T/300)^(0.0) * (875.0/(1-omega)),  # Reaction rate number 849 with type CR
         :k850 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 850 with type CR
         :k851 => 1.3e-17 * (T/300)^(0.0) * (875.0/(1-omega)),  # Reaction rate number 851 with type CR
         :k852 => 1.3e-17 * (T/300)^(0.0) * (2652.5/(1-omega)),  # Reaction rate number 852 with type CR
         :k853 => 1.3e-17 * (T/300)^(0.0) * (1642.0/(1-omega)),  # Reaction rate number 853 with type CR
         :k854 => 1.3e-17 * (T/300)^(0.0) * (1881.0/(1-omega)),  # Reaction rate number 854 with type CR
         :k855 => 1.3e-17 * (T/300)^(0.0) * (389.0/(1-omega)),  # Reaction rate number 855 with type CR
         :k856 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 856 with type CR
         :k857 => 1.3e-17 * (T/300)^(0.0) * (559.5/(1-omega)),  # Reaction rate number 857 with type CR
         :k858 => 1.3e-17 * (T/300)^(0.0) * (263.5/(1-omega)),  # Reaction rate number 858 with type CR
         :k859 => 1.3e-17 * (T/300)^(0.0) * (263.5/(1-omega)),  # Reaction rate number 859 with type CR
         :k860 => 1.3e-17 * (T/300)^(0.0) * (1122.5/(1-omega)),  # Reaction rate number 860 with type CR
         :k861 => 1.3e-17 * (T/300)^(0.0) * (2388.0/(1-omega)),  # Reaction rate number 861 with type CR
         :k862 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 862 with type CR
         :k863 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 863 with type CR
         :k864 => 1.3e-17 * (T/300)^(0.0) * (1000.0/(1-omega)),  # Reaction rate number 864 with type CR
         :k865 => 1.3e-17 * (T/300)^(0.0) * (559.5/(1-omega)),  # Reaction rate number 865 with type CR
         :k866 => 1.3e-17 * (T/300)^(0.0) * (857.0/(1-omega)),  # Reaction rate number 866 with type CR
         :k867 => 1.3e-17 * (T/300)^(0.0) * (717.5/(1-omega)),  # Reaction rate number 867 with type CR
         :k868 => 1.3e-17 * (T/300)^(0.0) * (1584.0/(1-omega)),  # Reaction rate number 868 with type CR
         :k869 => 1.3e-17 * (T/300)^(0.0) * (752.0/(1-omega)),  # Reaction rate number 869 with type CR
         :k870 => 1.3e-17 * (T/300)^(0.0) * (1169.5/(1-omega)),  # Reaction rate number 870 with type CR
         :k871 => 1.3e-17 * (T/300)^(0.0) * (365.0/(1-omega)),  # Reaction rate number 871 with type CR
         :k872 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 872 with type CR
         :k873 => 1.3e-17 * (T/300)^(0.0) * (5290.0/(1-omega)),  # Reaction rate number 873 with type CR
         :k874 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 874 with type CR
         :k875 => 1.3e-17 * (T/300)^(0.0) * (854.0/(1-omega)),  # Reaction rate number 875 with type CR
         :k876 => 1.3e-17 * (T/300)^(1.17) * (105.0/(1-omega)),  # Reaction rate number 876 with type CR
         :k877 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 877 with type CR
         :k878 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 878 with type CR
         :k879 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 879 with type CR
         :k880 => 1.3e-17 * (T/300)^(0.0) * (61.0/(1-omega)),  # Reaction rate number 880 with type CR
         :k881 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 881 with type CR
         :k882 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 882 with type CR
         :k883 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 883 with type CR
         :k884 => 1.3e-17 * (T/300)^(0.0) * (2500.0/(1-omega)),  # Reaction rate number 884 with type CR
         :k885 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 885 with type CR
         :k886 => 1.3e-17 * (T/300)^(0.0) * (1329.5/(1-omega)),  # Reaction rate number 886 with type CR
         :k887 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 887 with type CR
         :k888 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 888 with type CR
         :k889 => 1.3e-17 * (T/300)^(0.0) * (485.5/(1-omega)),  # Reaction rate number 889 with type CR
         :k890 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 890 with type CR
         :k891 => 1.3e-17 * (T/300)^(0.0) * (848.0/(1-omega)),  # Reaction rate number 891 with type CR
         :k892 => 1.3e-17 * (T/300)^(0.0) * (2577.0/(1-omega)),  # Reaction rate number 892 with type CR
         :k893 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 893 with type CR
         :k894 => 1.3e-17 * (T/300)^(0.0) * (0.2/(1-omega)),  # Reaction rate number 894 with type CR
         :k895 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 895 with type CR
         :k896 => 1.3e-17 * (T/300)^(0.0) * (863.5/(1-omega)),  # Reaction rate number 896 with type CR
         :k897 => 1.3e-17 * (T/300)^(0.0) * (865.0/(1-omega)),  # Reaction rate number 897 with type CR
         :k898 => 1.3e-17 * (T/300)^(0.0) * (559.5/(1-omega)),  # Reaction rate number 898 with type CR
         :k899 => 1.3e-17 * (T/300)^(0.0) * (865.0/(1-omega)),  # Reaction rate number 899 with type CR
         :k900 => 1.3e-17 * (T/300)^(0.0) * (875.0/(1-omega)),  # Reaction rate number 900 with type CR
         :k901 => 1.3e-17 * (T/300)^(0.0) * (875.0/(1-omega)),  # Reaction rate number 901 with type CR
         :k902 => 1.3e-17 * (T/300)^(0.0) * (875.0/(1-omega)),  # Reaction rate number 902 with type CR
         :k903 => 1.3e-17 * (T/300)^(0.0) * (875.0/(1-omega)),  # Reaction rate number 903 with type CR
         :k904 => 1.3e-17 * (T/300)^(0.0) * (1557.0/(1-omega)),  # Reaction rate number 904 with type CR
         :k905 => 1.3e-17 * (T/300)^(0.0) * (1500.0/(1-omega)),  # Reaction rate number 905 with type CR
         :k906 => 1.3e-17 * (T/300)^(0.0) * (210.5/(1-omega)),  # Reaction rate number 906 with type CR
         :k907 => 1.3e-17 * (T/300)^(0.0) * (584.5/(1-omega)),  # Reaction rate number 907 with type CR
         :k908 => 1.3e-17 * (T/300)^(0.0) * (1000.0/(1-omega)),  # Reaction rate number 908 with type CR
         :k909 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 909 with type CR
         :k910 => 1.3e-17 * (T/300)^(0.0) * (124.5/(1-omega)),  # Reaction rate number 910 with type CR
         :k911 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 911 with type CR
         :k912 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 912 with type CR
         :k913 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 913 with type CR
         :k914 => 1.3e-17 * (T/300)^(0.0) * (1370.0/(1-omega)),  # Reaction rate number 914 with type CR
         :k915 => 1.3e-17 * (T/300)^(0.0) * (124.0/(1-omega)),  # Reaction rate number 915 with type CR
         :k916 => 1.3e-17 * (T/300)^(0.0) * (1725.0/(1-omega)),  # Reaction rate number 916 with type CR
         :k917 => 1.3e-17 * (T/300)^(0.0) * (1500.0/(1-omega)),  # Reaction rate number 917 with type CR
         :k918 => 1.3e-17 * (T/300)^(0.0) * (1500.0/(1-omega)),  # Reaction rate number 918 with type CR
         :k919 => 1.3e-17 * (T/300)^(0.0) * (500.0/(1-omega)),  # Reaction rate number 919 with type CR
         :k920 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 920 with type CR
         :k921 => 1.3e-17 * (T/300)^(0.0) * (1500.0/(1-omega)),  # Reaction rate number 921 with type CR
         :k922 => 1.3e-17 * (T/300)^(0.0) * (1500.0/(1-omega)),  # Reaction rate number 922 with type CR
         :k923 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 923 with type CR
         :k924 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 924 with type CR
         :k925 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 925 with type CR
         :k926 => 1.3e-17 * (T/300)^(0.0) * (0.2/(1-omega)),  # Reaction rate number 926 with type CR
         :k927 => 1.3e-17 * (T/300)^(0.0) * (66.5/(1-omega)),  # Reaction rate number 927 with type CR
         :k928 => 1.3e-17 * (T/300)^(0.0) * (25.0/(1-omega)),  # Reaction rate number 928 with type CR
         :k929 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 929 with type CR
         :k930 => 1.3e-17 * (T/300)^(0.0) * (1.1/(1-omega)),  # Reaction rate number 930 with type CR
         :k931 => 1.3e-17 * (T/300)^(0.0) * (474.0/(1-omega)),  # Reaction rate number 931 with type CR
         :k932 => 1.3e-17 * (T/300)^(0.0) * (324.5/(1-omega)),  # Reaction rate number 932 with type CR
         :k933 => 1.3e-17 * (T/300)^(0.0) * (40.0/(1-omega)),  # Reaction rate number 933 with type CR
         :k934 => 1.3e-17 * (T/300)^(0.0) * (475.0/(1-omega)),  # Reaction rate number 934 with type CR
         :k935 => 1.3e-17 * (T/300)^(0.0) * (657.5/(1-omega)),  # Reaction rate number 935 with type CR
         :k936 => 1.3e-17 * (T/300)^(0.0) * (288.0/(1-omega)),  # Reaction rate number 936 with type CR
         :k937 => 1.3e-17 * (T/300)^(0.0) * (270.5/(1-omega)),  # Reaction rate number 937 with type CR
         :k938 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 938 with type CR
         :k939 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 939 with type CR
         :k940 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 940 with type CR
         :k941 => 1.3e-17 * (T/300)^(0.0) * (247.0/(1-omega)),  # Reaction rate number 941 with type CR
         :k942 => 1.3e-17 * (T/300)^(0.0) * (231.0/(1-omega)),  # Reaction rate number 942 with type CR
         :k943 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 943 with type CR
         :k944 => 1.3e-17 * (T/300)^(0.0) * (8.5/(1-omega)),  # Reaction rate number 944 with type CR
         :k945 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 945 with type CR
         :k946 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 946 with type CR
         :k947 => 1.3e-17 * (T/300)^(0.0) * (58.5/(1-omega)),  # Reaction rate number 947 with type CR
         :k948 => 1.3e-17 * (T/300)^(0.0) * (375.5/(1-omega)),  # Reaction rate number 948 with type CR
         :k949 => 1.3e-17 * (T/300)^(0.0) * (375.0/(1-omega)),  # Reaction rate number 949 with type CR
         :k950 => 1.3e-17 * (T/300)^(0.0) * (1.4/(1-omega)),  # Reaction rate number 950 with type CR
         :k951 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 951 with type CR
         :k952 => 1.3e-17 * (T/300)^(0.0) * (722.0/(1-omega)),  # Reaction rate number 952 with type CR
         :k953 => 1.3e-17 * (T/300)^(0.0) * (2680.0/(1-omega)),  # Reaction rate number 953 with type CR
         :k954 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 954 with type CR
         :k955 => 1.3e-17 * (T/300)^(0.0) * (254.5/(1-omega)),  # Reaction rate number 955 with type CR
         :k956 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 956 with type CR
         :k957 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 957 with type CR
         :k958 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 958 with type CR
         :k959 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 959 with type CR
         :k960 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 960 with type CR
         :k961 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 961 with type CR
         :k962 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 962 with type CR
         :k963 => 1.3e-17 * (T/300)^(0.0) * (480.0/(1-omega)),  # Reaction rate number 963 with type CR
         :k964 => 1.3e-17 * (T/300)^(0.0) * (922.0/(1-omega)),  # Reaction rate number 964 with type CR
         :k965 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 965 with type CR
         :k966 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 966 with type CR
         :k967 => 1.3e-17 * (T/300)^(0.0) * (2115.0/(1-omega)),  # Reaction rate number 967 with type CR
         :k968 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 968 with type CR
         :k969 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 969 with type CR
         :k970 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 970 with type CR
         :k971 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 971 with type CR
         :k972 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 972 with type CR
         :k973 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 973 with type CR
         :k974 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 974 with type CR
         :k975 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 975 with type CR
         :k976 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 976 with type CR
         :k977 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 977 with type CR
         :k978 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 978 with type CR
         :k979 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 979 with type CR
         :k980 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 980 with type CR
         :k981 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 981 with type CR
         :k982 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 982 with type CR
         :k983 => 1.3e-17 * (T/300)^(0.0) * (750.0/(1-omega)),  # Reaction rate number 983 with type CR
         :k984 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 984 with type CR
         :k985 => 1.3e-17 * (T/300)^(0.0) * (250.0/(1-omega)),  # Reaction rate number 985 with type CR
         :k986 => 1e-06 * (T/300)^(-0.3) * exp(-(0.0/T)),  # Reaction rate number 986 with type DR
         :k987 => 1e-06 * (T/300)^(-0.3) * exp(-(0.0/T)),  # Reaction rate number 987 with type DR
         :k988 => 1.53e-06 * (T/300)^(-0.3) * exp(-(0.0/T)),  # Reaction rate number 988 with type DR
         :k989 => 5.4e-08 * (T/300)^(-0.3) * exp(-(0.0/T)),  # Reaction rate number 989 with type DR
         :k990 => 3.92e-07 * (T/300)^(-0.3) * exp(-(0.0/T)),  # Reaction rate number 990 with type DR
         :k991 => 1.35e-06 * (T/300)^(-0.3) * exp(-(0.0/T)),  # Reaction rate number 991 with type DR
         :k992 => 1e-06 * (T/300)^(-0.3) * exp(-(0.0/T)),  # Reaction rate number 992 with type DR
         :k993 => 1e-06 * (T/300)^(-0.3) * exp(-(0.0/T)),  # Reaction rate number 993 with type DR
         :k994 => 8.39e-08 * (T/300)^(-0.3) * exp(-(0.0/T)),  # Reaction rate number 994 with type DR
         :k995 => 2.66e-08 * (T/300)^(-0.3) * exp(-(0.0/T)),  # Reaction rate number 995 with type DR
         :k996 => 1.2e-06 * (T/300)^(-0.3) * exp(-(0.0/T)),  # Reaction rate number 996 with type DR
         :k997 => 6.76e-07 * (T/300)^(-0.3) * exp(-(0.0/T)),  # Reaction rate number 997 with type DR
         :k998 => 3e-07 * (T/300)^(-0.5) * exp(-(0.0/T)),  # Reaction rate number 998 with type DR
         :k999 => 1.16e-07 * (T/300)^(-0.76) * exp(-(0.0/T)),  # Reaction rate number 999 with type DR
         :k1000 => 1.53e-07 * (T/300)^(-0.76) * exp(-(0.0/T)) # Reaction rate number 1000 with type DR
         )


### Tempurature range flags ###
T_low_bounds = [10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 88 10 10 10 10 10 10 10 10 1340 2803 1763 20 1743 1696 3400 1340 1833 1763 1743 1696 200 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 205 10 10 10 10 200 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 2000 10 10 10 10 10 10 10 10 10 10 10 10 10 200 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 ]
T_upp_bounds = [41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 41000 41000 41000 41000 41000 41000 41000 2500 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 100 41000 41000 41000 41000 41000 41000 41000 41000 41000 300 300 300 300 300 300 41000 300 300 300 300 41000 300 300 300 300 300 300 300 300 300 41000 41000 41000 41000 300 300 300 300 300 300 300 300 300 300 300 300 41000 41000 41000 494 41000 41000 300 300 300 300 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 700 41000 41000 41000 300 41000 300 41000 300 300 300 41000 41000 41000 41000 41000 300 300 41000 41000 41000 300 41000 300 41000 41000 41000 41000 41000 41000 300 300 41000 300 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 300 41000 300 41000 41000 300 300 300 300 300 300 300 300 300 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 300 41000 300 300 300 41000 300 41000 41000 41000 300 41000 41000 41000 41000 41000 41000 41000 300 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 300 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 300 41000 300 41000 300 41000 300 41000 41000 41000 41000 300 41000 300 41000 41000 300 300 41000 41000 41000 300 41000 300 300 300 300 300 300 300 300 300 41000 41000 41000 300 41000 300 300 41000 300 41000 300 300 41000 41000 41000 41000 41000 41000 300 41000 300 300 41000 300 41000 41000 41000 41000 300 300 565 41000 300 41000 300 520 300 41000 300 41000 41000 41000 41000 41000 41000 300 41000 300 41000 300 41000 300 41000 41000 41000 300 300 41000 300 41000 41000 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 300 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 300 41000 41000 41000 41000 41000 41000 41000 41000 300 41000 41000 41000 300 41000 41000 300 300 300 41000 41000 41000 300 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 300 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 300 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 300 41000 41000 300 41000 41000 10000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 700 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 4999 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 300 300 41000 41000 41000 41000 41000 41000 41000 41000 41000 41000 300 300 300 ]
for i = 1:length(T_low_bounds)
    if T < T_low_bounds[i]
        print("\nTempurature below the tempurature range [", T_low_bounds[i], ", ", T_upp_bounds[i], "] for reaction number ", i ) 
    elseif T > T_upp_bounds[i]
        print("\nTempurature higher than the tempurature range [", T_low_bounds[i], ", ", T_upp_bounds[i], "] for reaction number ", i ) 
    end
end
print("\n")
print("\n")


### Network Admin things ###
@variables t 
@species C⁻(t) C(t) C2(t) e(t) CH2(t) C2H2(t) CH(t) C2H(t) CO2(t) CO(t) H2O(t) H2CO(t) N(t) CN(t) NH(t) HCN(t) O2(t) O(t) OH(t) HCO(t) C2⁻(t) C4(t) C3(t) C5(t) C6(t) C7(t) C8(t) C9(t) C10(t) C3⁻(t) C4⁻(t) C5⁻(t) C6⁻(t) C7⁻(t) C8⁻(t) C10⁻(t) C2H⁻(t) C3H(t) C3H⁻(t) C4H(t) C3N⁻(t) C4N(t) C4H⁻(t) C5H(t) C5H⁻(t) C6H(t) C6H⁻(t) C7H(t) C7H⁻(t) C8H(t) C8H⁻(t) C9H(t) C9⁻(t) C9H⁻(t) C10H(t) CH⁻(t) O⁻(t) OH⁻(t) S⁻(t) CS(t) CH3(t) CN⁻(t) CH3CN(t) CH3OH(t) HCO⁺(t) OCS(t) H⁻(t) CH4(t) H(t) H2(t) NH2(t) NH3(t) C10H⁻(t) C10H2(t) H2CCC(t) HC3N(t) HC4H(t) C5H2(t) C5N⁻(t) HC5N(t) C6H2(t) C7H2(t) C8H2(t) C9H2(t) HS(t) C2N(t) C3N(t) C5N(t) C7N(t) HC7N(t) C9N(t) HC9N(t) NO(t) NS(t) NO2(t) SO2(t) C2O(t) SO(t) H⁺(t) HNC(t) HOC⁺(t) O2⁻(t) C⁺(t) C10H⁺(t) C2H4(t) C2H4⁺(t) C2H5(t) C2H5⁺(t) C2H5OH(t) C2H5OH⁺(t) C2O⁺(t) C2S(t) C2S⁺(t) C3O(t) C3O⁺(t) C3S(t) C3S⁺(t) C4H3(t) C4H3⁺(t) C4S(t) C4S⁺(t) C6H6(t) C6H6⁺(t) CCP(t) CCP⁺(t) CCl(t) CCl⁺(t) CH2⁺(t) CH2CCH2(t) C3H4⁺(t) CH2CCH(t) CH2CCH⁺(t) CH2CN(t) CH2CN⁺(t) CH2CO(t) CH2CO⁺(t) CH3CCH(t) CH3CHCH2(t) C3H6⁺(t) CH3CHO(t) CH3CHO⁺(t) CH3COCH3(t) CH3COCH3⁺(t) CH3OCH3(t) CH3OCH3⁺(t) CH⁺(t) CP(t) CP⁺(t) ClO(t) ClO⁺(t) Fe(t) Fe⁺(t) H2CO⁺(t) H2S(t) H2S⁺(t) H2SiO(t) H2SiO⁺(t) C4H2⁺(t) HCOOCH3(t) COOCH4⁺(t) HCP(t) HCP⁺(t) HPO(t) HPO⁺(t) Mg(t) Mg⁺(t) NCCN(t) C2N⁺(t) CNC⁺(t) NH3⁺(t) NO⁺(t) NS⁺(t) Na(t) Na⁺(t) OCS⁺(t) P(t) P⁺(t) PH(t) PH⁺(t) PO(t) PO⁺(t) SO⁺(t) Si(t) Si⁺(t) SiC2(t) SiC2⁺(t) SiC2H(t) SiC2H⁺(t) SiC3(t) SiC3⁺(t) SiC(t) SiC⁺(t) SiCH2(t) SiCH2⁺(t) SiCH3(t) SiCH3⁺(t) SiH2(t) SiH2⁺(t) SiH3(t) SiH3⁺(t) SiN(t) SiN⁺(t) SiS(t) SiS⁺(t) C2⁺(t) S(t) S⁺(t) CN⁺(t) CO⁺(t) N2⁺(t) N2(t) O2⁺(t) C2H⁺(t) C2H2⁺(t) C2H3(t) C2H3⁺(t) C5H2⁺(t) C6H2⁺(t) C7H2⁺(t) C3H3⁺(t) CO2⁺(t) HC3N⁺(t) HCN⁺(t) C2N2⁺(t) C3⁺(t) C5⁺(t) CH3CH3⁺(t) CH3CH3(t) PN⁺(t) PN(t) H2O⁺(t) NH2⁺(t) O⁺(t) OH⁺(t) CH3⁺(t) CH4⁺(t) CH3OH⁺(t) N⁺(t) SO2⁺(t) CS⁺(t) Cl⁺(t) Cl(t) C10⁺(t) C10H2⁺(t) C3H2(t) C3H2⁺(t) C3H⁺(t) C4⁺(t) C4H⁺(t) C4P(t) C4P⁺(t) C5H⁺(t) C6⁺(t) C6H⁺(t) C7⁺(t) C7H⁺(t) C8⁺(t) C8H2⁺(t) C8H⁺(t) C9⁺(t) C9H2⁺(t) C9H⁺(t) CH3C4H(t) CH3C4H⁺(t) CH3C6H(t) C7H4⁺(t) CH3CN⁺(t) H2CS(t) H2CS⁺(t) H2S2(t) H2S2⁺(t) HC2P(t) HC2P⁺(t) HC5N⁺(t) HC7N⁺(t) HC9N⁺(t) HCSi(t) HCSi⁺(t) HCl(t) HCl⁺(t) HNSi(t) HNSi⁺(t) HS2(t) HS2⁺(t) HS⁺(t) N2O(t) N2O⁺(t) NH⁺(t) PH2(t) PH2⁺(t) S2(t) S2⁺(t) SiC2H2(t) SiC2H2⁺(t) SiC3H(t) SiC3H⁺(t) SiC4(t) SiC4⁺(t) SiH4(t) SiH4⁺(t) SiH(t) SiH⁺(t) SiNC(t) SiNC⁺(t) SiO(t) SiO⁺(t) H2⁺(t) F⁺(t) F(t) He⁺(t) He(t) C3N⁺(t) HCNO(t) HCNO⁺(t) HNCO(t) HNCO⁺(t) HONC(t) HONC⁺(t) HNO⁺(t) HNO(t) HCOOH(t) HCOOH⁺(t) NO2⁺(t) C11(t) C2H5CN(t) CH2CHCNH⁺(t) C3P(t) CH2CHCCH(t) CH2CHCHCH2(t) CH2CHCN(t) CH2NH(t) CH2PH(t) CH3C3N(t) CH3C5N(t) CH3C7N(t) CNO(t) H2CN(t) H2O2(t) HCS(t) HCS⁺(t) HF(t) HNC3(t) HOCN(t) NH2CN(t) O2H(t) OCN(t) SiO2(t) C10H3⁺(t) C11⁺(t) 
@parameters T Av n_H cr_ion_rate omega k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14 k15 k16 k17 k18 k19 k20 k21 k22 k23 k24 k25 k26 k27 k28 k29 k30 k31 k32 k33 k34 k35 k36 k37 k38 k39 k40 k41 k42 k43 k44 k45 k46 k47 k48 k49 k50 k51 k52 k53 k54 k55 k56 k57 k58 k59 k60 k61 k62 k63 k64 k65 k66 k67 k68 k69 k70 k71 k72 k73 k74 k75 k76 k77 k78 k79 k80 k81 k82 k83 k84 k85 k86 k87 k88 k89 k90 k91 k92 k93 k94 k95 k96 k97 k98 k99 k100 k101 k102 k103 k104 k105 k106 k107 k108 k109 k110 k111 k112 k113 k114 k115 k116 k117 k118 k119 k120 k121 k122 k123 k124 k125 k126 k127 k128 k129 k130 k131 k132 k133 k134 k135 k136 k137 k138 k139 k140 k141 k142 k143 k144 k145 k146 k147 k148 k149 k150 k151 k152 k153 k154 k155 k156 k157 k158 k159 k160 k161 k162 k163 k164 k165 k166 k167 k168 k169 k170 k171 k172 k173 k174 k175 k176 k177 k178 k179 k180 k181 k182 k183 k184 k185 k186 k187 k188 k189 k190 k191 k192 k193 k194 k195 k196 k197 k198 k199 k200 k201 k202 k203 k204 k205 k206 k207 k208 k209 k210 k211 k212 k213 k214 k215 k216 k217 k218 k219 k220 k221 k222 k223 k224 k225 k226 k227 k228 k229 k230 k231 k232 k233 k234 k235 k236 k237 k238 k239 k240 k241 k242 k243 k244 k245 k246 k247 k248 k249 k250 k251 k252 k253 k254 k255 k256 k257 k258 k259 k260 k261 k262 k263 k264 k265 k266 k267 k268 k269 k270 k271 k272 k273 k274 k275 k276 k277 k278 k279 k280 k281 k282 k283 k284 k285 k286 k287 k288 k289 k290 k291 k292 k293 k294 k295 k296 k297 k298 k299 k300 k301 k302 k303 k304 k305 k306 k307 k308 k309 k310 k311 k312 k313 k314 k315 k316 k317 k318 k319 k320 k321 k322 k323 k324 k325 k326 k327 k328 k329 k330 k331 k332 k333 k334 k335 k336 k337 k338 k339 k340 k341 k342 k343 k344 k345 k346 k347 k348 k349 k350 k351 k352 k353 k354 k355 k356 k357 k358 k359 k360 k361 k362 k363 k364 k365 k366 k367 k368 k369 k370 k371 k372 k373 k374 k375 k376 k377 k378 k379 k380 k381 k382 k383 k384 k385 k386 k387 k388 k389 k390 k391 k392 k393 k394 k395 k396 k397 k398 k399 k400 k401 k402 k403 k404 k405 k406 k407 k408 k409 k410 k411 k412 k413 k414 k415 k416 k417 k418 k419 k420 k421 k422 k423 k424 k425 k426 k427 k428 k429 k430 k431 k432 k433 k434 k435 k436 k437 k438 k439 k440 k441 k442 k443 k444 k445 k446 k447 k448 k449 k450 k451 k452 k453 k454 k455 k456 k457 k458 k459 k460 k461 k462 k463 k464 k465 k466 k467 k468 k469 k470 k471 k472 k473 k474 k475 k476 k477 k478 k479 k480 k481 k482 k483 k484 k485 k486 k487 k488 k489 k490 k491 k492 k493 k494 k495 k496 k497 k498 k499 k500 k501 k502 k503 k504 k505 k506 k507 k508 k509 k510 k511 k512 k513 k514 k515 k516 k517 k518 k519 k520 k521 k522 k523 k524 k525 k526 k527 k528 k529 k530 k531 k532 k533 k534 k535 k536 k537 k538 k539 k540 k541 k542 k543 k544 k545 k546 k547 k548 k549 k550 k551 k552 k553 k554 k555 k556 k557 k558 k559 k560 k561 k562 k563 k564 k565 k566 k567 k568 k569 k570 k571 k572 k573 k574 k575 k576 k577 k578 k579 k580 k581 k582 k583 k584 k585 k586 k587 k588 k589 k590 k591 k592 k593 k594 k595 k596 k597 k598 k599 k600 k601 k602 k603 k604 k605 k606 k607 k608 k609 k610 k611 k612 k613 k614 k615 k616 k617 k618 k619 k620 k621 k622 k623 k624 k625 k626 k627 k628 k629 k630 k631 k632 k633 k634 k635 k636 k637 k638 k639 k640 k641 k642 k643 k644 k645 k646 k647 k648 k649 k650 k651 k652 k653 k654 k655 k656 k657 k658 k659 k660 k661 k662 k663 k664 k665 k666 k667 k668 k669 k670 k671 k672 k673 k674 k675 k676 k677 k678 k679 k680 k681 k682 k683 k684 k685 k686 k687 k688 k689 k690 k691 k692 k693 k694 k695 k696 k697 k698 k699 k700 k701 k702 k703 k704 k705 k706 k707 k708 k709 k710 k711 k712 k713 k714 k715 k716 k717 k718 k719 k720 k721 k722 k723 k724 k725 k726 k727 k728 k729 k730 k731 k732 k733 k734 k735 k736 k737 k738 k739 k740 k741 k742 k743 k744 k745 k746 k747 k748 k749 k750 k751 k752 k753 k754 k755 k756 k757 k758 k759 k760 k761 k762 k763 k764 k765 k766 k767 k768 k769 k770 k771 k772 k773 k774 k775 k776 k777 k778 k779 k780 k781 k782 k783 k784 k785 k786 k787 k788 k789 k790 k791 k792 k793 k794 k795 k796 k797 k798 k799 k800 k801 k802 k803 k804 k805 k806 k807 k808 k809 k810 k811 k812 k813 k814 k815 k816 k817 k818 k819 k820 k821 k822 k823 k824 k825 k826 k827 k828 k829 k830 k831 k832 k833 k834 k835 k836 k837 k838 k839 k840 k841 k842 k843 k844 k845 k846 k847 k848 k849 k850 k851 k852 k853 k854 k855 k856 k857 k858 k859 k860 k861 k862 k863 k864 k865 k866 k867 k868 k869 k870 k871 k872 k873 k874 k875 k876 k877 k878 k879 k880 k881 k882 k883 k884 k885 k886 k887 k888 k889 k890 k891 k892 k893 k894 k895 k896 k897 k898 k899 k900 k901 k902 k903 k904 k905 k906 k907 k908 k909 k910 k911 k912 k913 k914 k915 k916 k917 k918 k919 k920 k921 k922 k923 k924 k925 k926 k927 k928 k929 k930 k931 k932 k933 k934 k935 k936 k937 k938 k939 k940 k941 k942 k943 k944 k945 k946 k947 k948 k949 k950 k951 k952 k953 k954 k955 k956 k957 k958 k959 k960 k961 k962 k963 k964 k965 k966 k967 k968 k969 k970 k971 k972 k973 k974 k975 k976 k977 k978 k979 k980 k981 k982 k983 k984 k985 k986 k987 k988 k989 k990 k991 k992 k993 k994 k995 k996 k997 k998 k999 k1000 


### Reaction Equations ###
reaction_equations = [
    (@reaction n_H * k1, C⁻ + C --> C2 + e ), 
    (@reaction n_H * k2, C⁻ + CH2 --> C2H2 + e ), 
    (@reaction n_H * k3, C⁻ + CH --> C2H + e ), 
    (@reaction n_H * k4, C⁻ + CO2 --> CO + CO + e ), 
    (@reaction n_H * k5, C⁻ + H2O --> H2CO + e ), 
    (@reaction n_H * k6, C⁻ + N --> CN + e ), 
    (@reaction n_H * k7, C⁻ + NH --> HCN + e ), 
    (@reaction n_H * k8, C⁻ + O2 --> CO2 + e ), 
    (@reaction n_H * k9, C⁻ + O --> CO + e ), 
    (@reaction n_H * k10, C⁻ + OH --> HCO + e ), 
    (@reaction n_H * k11, C2⁻ + C2 --> C4 + e ), 
    (@reaction n_H * k12, C2⁻ + C3 --> C5 + e ), 
    (@reaction n_H * k13, C2⁻ + C4 --> C6 + e ), 
    (@reaction n_H * k14, C2⁻ + C5 --> C7 + e ), 
    (@reaction n_H * k15, C2⁻ + C6 --> C8 + e ), 
    (@reaction n_H * k16, C2⁻ + C7 --> C9 + e ), 
    (@reaction n_H * k17, C2⁻ + C8 --> C10 + e ), 
    (@reaction n_H * k18, C3⁻ + C3 --> C6 + e ), 
    (@reaction n_H * k19, C3⁻ + C4 --> C7 + e ), 
    (@reaction n_H * k20, C3⁻ + C5 --> C8 + e ), 
    (@reaction n_H * k21, C3⁻ + C6 --> C9 + e ), 
    (@reaction n_H * k22, C3⁻ + C7 --> C10 + e ), 
    (@reaction n_H * k23, C4⁻ + C2 --> C6 + e ), 
    (@reaction n_H * k24, C4⁻ + C3 --> C7 + e ), 
    (@reaction n_H * k25, C4⁻ + C4 --> C8 + e ), 
    (@reaction n_H * k26, C4⁻ + C5 --> C9 + e ), 
    (@reaction n_H * k27, C4⁻ + C6 --> C10 + e ), 
    (@reaction n_H * k28, C5⁻ + C2 --> C7 + e ), 
    (@reaction n_H * k29, C5⁻ + C3 --> C8 + e ), 
    (@reaction n_H * k30, C5⁻ + C4 --> C9 + e ), 
    (@reaction n_H * k31, C5⁻ + C5 --> C10 + e ), 
    (@reaction n_H * k32, C6⁻ + C2 --> C8 + e ), 
    (@reaction n_H * k33, C6⁻ + C3 --> C9 + e ), 
    (@reaction n_H * k34, C6⁻ + C4 --> C10 + e ), 
    (@reaction n_H * k35, C7⁻ + C2 --> C9 + e ), 
    (@reaction n_H * k36, C7⁻ + C3 --> C10 + e ), 
    (@reaction n_H * k37, C8⁻ + C2 --> C10 + e ), 
    (@reaction n_H * k38, C + C10⁻ --> C5 + C6 + e ), 
    (@reaction n_H * k39, C + C2⁻ --> C3 + e ), 
    (@reaction n_H * k40, C + C2H⁻ --> C3H + e ), 
    (@reaction n_H * k41, C + C3⁻ --> C4 + e ), 
    (@reaction n_H * k42, C + C3H⁻ --> C4H + e ), 
    (@reaction n_H * k43, C + C3N⁻ --> C4N + e ), 
    (@reaction n_H * k44, C + C4⁻ --> C5 + e ), 
    (@reaction n_H * k45, C + C4H⁻ --> C5H + e ), 
    (@reaction n_H * k46, C + C5⁻ --> C6 + e ), 
    (@reaction n_H * k47, C + C5H⁻ --> C6H + e ), 
    (@reaction n_H * k48, C + C6⁻ --> C7 + e ), 
    (@reaction n_H * k49, C + C6H⁻ --> C7H + e ), 
    (@reaction n_H * k50, C + C7⁻ --> C8 + e ), 
    (@reaction n_H * k51, C + C7H⁻ --> C8H + e ), 
    (@reaction n_H * k52, C + C8⁻ --> C9 + e ), 
    (@reaction n_H * k53, C + C8H⁻ --> C9H + e ), 
    (@reaction n_H * k54, C + C9⁻ --> C10 + e ), 
    (@reaction n_H * k55, C + C9H⁻ --> C10H + e ), 
    (@reaction n_H * k56, C + CH⁻ --> C2H + e ), 
    (@reaction n_H * k57, C + O⁻ --> CO + e ), 
    (@reaction n_H * k58, C + OH⁻ --> HCO + e ), 
    (@reaction n_H * k59, C + S⁻ --> CS + e ), 
    (@reaction n_H * k60, CH2 + O⁻ --> H2CO + e ), 
    (@reaction n_H * k61, CH3 + CN⁻ --> CH3CN + e ), 
    (@reaction n_H * k62, CH3 + OH⁻ --> CH3OH + e ), 
    (@reaction n_H * k63, CH + O⁻ --> HCO + e ), 
    (@reaction n_H * k64, CH + O --> HCO⁺ + e ), 
    (@reaction n_H * k65, CH + OH⁻ --> H2CO + e ), 
    (@reaction n_H * k66, CO + S⁻ --> OCS + e ), 
    (@reaction n_H * k67, H⁻ + C2 --> C2H + e ), 
    (@reaction n_H * k68, H⁻ + C2H --> C2H2 + e ), 
    (@reaction n_H * k69, H⁻ + C --> CH + e ), 
    (@reaction n_H * k70, H⁻ + CH2 --> CH3 + e ), 
    (@reaction n_H * k71, H⁻ + CH3 --> CH4 + e ), 
    (@reaction n_H * k72, H⁻ + CH --> CH2 + e ), 
    (@reaction n_H * k73, H⁻ + CN --> HCN + e ), 
    (@reaction n_H * k74, H⁻ + CO --> HCO + e ), 
    (@reaction n_H * k75, H⁻ + H --> H2 + e ), 
    (@reaction n_H * k76, H⁻ + HCO --> H2CO + e ), 
    (@reaction n_H * k77, H⁻ + N --> NH + e ), 
    (@reaction n_H * k78, H⁻ + NH2 --> NH3 + e ), 
    (@reaction n_H * k79, H⁻ + NH --> NH2 + e ), 
    (@reaction n_H * k80, H⁻ + O --> OH + e ), 
    (@reaction n_H * k81, H⁻ + OH --> H2O + e ), 
    (@reaction n_H * k82, H2 + C⁻ --> CH2 + e ), 
    (@reaction n_H * k83, H2 + O⁻ --> H2O + e ), 
    (@reaction n_H * k84, H + C⁻ --> CH + e ), 
    (@reaction n_H * k85, H + C10⁻ --> C10H + e ), 
    (@reaction n_H * k86, H + C10H⁻ --> C10H2 + e ), 
    (@reaction n_H * k87, H + C2⁻ --> C2H + e ), 
    (@reaction n_H * k88, H + C2H⁻ --> C2H2 + e ), 
    (@reaction n_H * k89, H + C3⁻ --> C3H + e ), 
    (@reaction n_H * k90, H + C3H⁻ --> H2CCC + e ), 
    (@reaction n_H * k91, H + C3N⁻ --> HC3N + e ), 
    (@reaction n_H * k92, H + C4⁻ --> C4H + e ), 
    (@reaction n_H * k93, H + C4H⁻ --> HC4H + e ), 
    (@reaction n_H * k94, H + C5⁻ --> C5H + e ), 
    (@reaction n_H * k95, H + C5H⁻ --> C5H2 + e ), 
    (@reaction n_H * k96, H + C5N⁻ --> HC5N + e ), 
    (@reaction n_H * k97, H + C6⁻ --> C6H + e ), 
    (@reaction n_H * k98, H + C6H⁻ --> C6H2 + e ), 
    (@reaction n_H * k99, H + C7⁻ --> C7H + e ), 
    (@reaction n_H * k100, H + C7H⁻ --> C7H2 + e ), 
    (@reaction n_H * k101, H + C8⁻ --> C8H + e ), 
    (@reaction n_H * k102, H + C8H⁻ --> C8H2 + e ), 
    (@reaction n_H * k103, H + C9⁻ --> C9H + e ), 
    (@reaction n_H * k104, H + C9H⁻ --> C9H2 + e ), 
    (@reaction n_H * k105, H + CH⁻ --> CH2 + e ), 
    (@reaction n_H * k106, H + CN⁻ --> HCN + e ), 
    (@reaction n_H * k107, H + O⁻ --> OH + e ), 
    (@reaction n_H * k108, H + OH⁻ --> H2O + e ), 
    (@reaction n_H * k109, H + S⁻ --> HS + e ), 
    (@reaction n_H * k110, N + C2⁻ --> C2N + e ), 
    (@reaction n_H * k111, N + C2H⁻ --> C2N + H + e ), 
    (@reaction n_H * k112, N + C3⁻ --> C3N + e ), 
    (@reaction n_H * k113, N + C3H⁻ --> HC3N + e ), 
    (@reaction n_H * k114, N + C4⁻ --> C4N + e ), 
    (@reaction n_H * k115, N + C4H⁻ --> CN + C3H + e ), 
    (@reaction n_H * k116, N + C5⁻ --> C5N + e ), 
    (@reaction n_H * k117, N + C5H⁻ --> HC5N + e ), 
    (@reaction n_H * k118, N + C7⁻ --> C7N + e ), 
    (@reaction n_H * k119, N + C7H⁻ --> HC7N + e ), 
    (@reaction n_H * k120, N + C9⁻ --> C9N + e ), 
    (@reaction n_H * k121, N + C9H⁻ --> HC9N + e ), 
    (@reaction n_H * k122, N + O⁻ --> NO + e ), 
    (@reaction n_H * k123, N + S⁻ --> NS + e ), 
    (@reaction n_H * k124, O⁻ + CO --> CO2 + e ), 
    (@reaction n_H * k125, O⁻ + NO --> NO2 + e ), 
    (@reaction n_H * k126, O⁻ + O --> O2 + e ), 
    (@reaction n_H * k127, O2 + S⁻ --> SO2 + e ), 
    (@reaction n_H * k128, O + C2⁻ --> C2O + e ), 
    (@reaction n_H * k129, O + C2H⁻ --> C2O + H + e ), 
    (@reaction n_H * k130, O + C3N⁻ --> CO + C2N + e ), 
    (@reaction n_H * k131, O + C5N⁻ --> CO + C4N + e ), 
    (@reaction n_H * k132, O + S⁻ --> SO + e ), 
    (@reaction n_H * k133, H⁺ + HNC --> HCN + H⁺ ), 
    (@reaction n_H * k134, H2 + CH --> C + H2 + H ), 
    (@reaction n_H * k135, H2 + H2 --> H2 + H + H ), 
    (@reaction n_H * k136, H2 + H2O --> OH + H2 + H ), 
    (@reaction n_H * k137, H2 + HOC⁺ --> HCO⁺ + H2 ), 
    (@reaction n_H * k138, H2 + O2 --> O + O + H2 ), 
    (@reaction n_H * k139, H2 + OH --> O + H2 + H ), 
    (@reaction n_H * k140, H2 + e --> H + H + e ), 
    (@reaction n_H * k141, H + CH --> C + H + H ), 
    (@reaction n_H * k142, H + H2 --> H + H + H ), 
    (@reaction n_H * k143, H + H2O --> OH + H + H ), 
    (@reaction n_H * k144, H + O2 --> O + O + H ), 
    (@reaction n_H * k145, H + OH --> O + H + H ), 
    (@reaction n_H * k146, O2⁻ + O2 --> O2 + O2 + e ), 
    (@reaction n_H * k147, C⁺ + C10H --> C10H⁺ + C ), 
    (@reaction n_H * k148, C⁺ + C2H4 --> C2H4⁺ + C ), 
    (@reaction n_H * k149, C⁺ + C2H5 --> C2H5⁺ + C ), 
    (@reaction n_H * k150, C⁺ + C2H5OH --> C2H5OH⁺ + C ), 
    (@reaction n_H * k151, C⁺ + C2O --> C2O⁺ + C ), 
    (@reaction n_H * k152, C⁺ + C2S --> C2S⁺ + C ), 
    (@reaction n_H * k153, C⁺ + C3O --> C3O⁺ + C ), 
    (@reaction n_H * k154, C⁺ + C3S --> C3S⁺ + C ), 
    (@reaction n_H * k155, C⁺ + C4H3 --> C4H3⁺ + C ), 
    (@reaction n_H * k156, C⁺ + C4S --> C4S⁺ + C ), 
    (@reaction n_H * k157, C⁺ + C6H6 --> C6H6⁺ + C ), 
    (@reaction n_H * k158, C⁺ + CCP --> CCP⁺ + C ), 
    (@reaction n_H * k159, C⁺ + CCl --> CCl⁺ + C ), 
    (@reaction n_H * k160, C⁺ + CH2 --> CH2⁺ + C ), 
    (@reaction n_H * k161, C⁺ + CH2CCH2 --> C3H4⁺ + C ), 
    (@reaction n_H * k162, C⁺ + CH2CCH --> CH2CCH⁺ + C ), 
    (@reaction n_H * k163, C⁺ + CH2CN --> CH2CN⁺ + C ), 
    (@reaction n_H * k164, C⁺ + CH2CO --> CH2CO⁺ + C ), 
    (@reaction n_H * k165, C⁺ + CH3CCH --> C3H4⁺ + C ), 
    (@reaction n_H * k166, C⁺ + CH3CHCH2 --> C3H6⁺ + C ), 
    (@reaction n_H * k167, C⁺ + CH3CHO --> CH3CHO⁺ + C ), 
    (@reaction n_H * k168, C⁺ + CH3COCH3 --> CH3COCH3⁺ + C ), 
    (@reaction n_H * k169, C⁺ + CH3OCH3 --> CH3OCH3⁺ + C ), 
    (@reaction n_H * k170, C⁺ + CH --> CH⁺ + C ), 
    (@reaction n_H * k171, C⁺ + CP --> CP⁺ + C ), 
    (@reaction n_H * k172, C⁺ + ClO --> ClO⁺ + C ), 
    (@reaction n_H * k173, C⁺ + Fe --> Fe⁺ + C ), 
    (@reaction n_H * k174, C⁺ + H2CO --> H2CO⁺ + C ), 
    (@reaction n_H * k175, C⁺ + H2S --> H2S⁺ + C ), 
    (@reaction n_H * k176, C⁺ + H2SiO --> H2SiO⁺ + C ), 
    (@reaction n_H * k177, C⁺ + HC4H --> C4H2⁺ + C ), 
    (@reaction n_H * k178, C⁺ + HCO --> HCO⁺ + C ), 
    (@reaction n_H * k179, C⁺ + HCOOCH3 --> COOCH4⁺ + C ), 
    (@reaction n_H * k180, C⁺ + HCP --> HCP⁺ + C ), 
    (@reaction n_H * k181, C⁺ + HPO --> HPO⁺ + C ), 
    (@reaction n_H * k182, C⁺ + Mg --> Mg⁺ + C ), 
    (@reaction n_H * k183, C⁺ + NCCN --> C2N⁺ + CN ), 
    (@reaction n_H * k184, C⁺ + NCCN --> CNC⁺ + CN ), 
    (@reaction n_H * k185, C⁺ + NH3 --> NH3⁺ + C ), 
    (@reaction n_H * k186, C⁺ + NO --> NO⁺ + C ), 
    (@reaction n_H * k187, C⁺ + NS --> NS⁺ + C ), 
    (@reaction n_H * k188, C⁺ + Na --> Na⁺ + C ), 
    (@reaction n_H * k189, C⁺ + OCS --> OCS⁺ + C ), 
    (@reaction n_H * k190, C⁺ + P --> P⁺ + C ), 
    (@reaction n_H * k191, C⁺ + PH --> PH⁺ + C ), 
    (@reaction n_H * k192, C⁺ + PO --> PO⁺ + C ), 
    (@reaction n_H * k193, C⁺ + SO --> SO⁺ + C ), 
    (@reaction n_H * k194, C⁺ + Si --> Si⁺ + C ), 
    (@reaction n_H * k195, C⁺ + SiC2 --> SiC2⁺ + C ), 
    (@reaction n_H * k196, C⁺ + SiC2H --> SiC2H⁺ + C ), 
    (@reaction n_H * k197, C⁺ + SiC3 --> SiC3⁺ + C ), 
    (@reaction n_H * k198, C⁺ + SiC --> SiC⁺ + C ), 
    (@reaction n_H * k199, C⁺ + SiCH2 --> SiCH2⁺ + C ), 
    (@reaction n_H * k200, C⁺ + SiCH3 --> SiCH3⁺ + C ), 
    (@reaction n_H * k201, C⁺ + SiH2 --> SiH2⁺ + C ), 
    (@reaction n_H * k202, C⁺ + SiH3 --> SiH3⁺ + C ), 
    (@reaction n_H * k203, C⁺ + SiN --> SiN⁺ + C ), 
    (@reaction n_H * k204, C⁺ + SiS --> SiS⁺ + C ), 
    (@reaction n_H * k205, C2⁺ + HCO --> HCO⁺ + C2 ), 
    (@reaction n_H * k206, C2⁺ + NO --> NO⁺ + C2 ), 
    (@reaction n_H * k207, C2⁺ + S --> S⁺ + C2 ), 
    (@reaction n_H * k208, C2 + CN⁺ --> CN + C2⁺ ), 
    (@reaction n_H * k209, C2 + CO⁺ --> CO + C2⁺ ), 
    (@reaction n_H * k210, C2 + N2⁺ --> N2 + C2⁺ ), 
    (@reaction n_H * k211, C2 + O2⁺ --> O2 + C2⁺ ), 
    (@reaction n_H * k212, C2H⁺ + NO --> NO⁺ + C2H ), 
    (@reaction n_H * k213, C2H⁺ + S --> S⁺ + C2H ), 
    (@reaction n_H * k214, C2H2⁺ + C2H3 --> C2H3⁺ + C2H2 ), 
    (@reaction n_H * k215, C2H2⁺ + C2H4 --> C2H4⁺ + C2H2 ), 
    (@reaction n_H * k216, C2H2⁺ + C5H2 --> C5H2⁺ + C2H2 ), 
    (@reaction n_H * k217, C2H2⁺ + C6H2 --> C6H2⁺ + C2H2 ), 
    (@reaction n_H * k218, C2H2⁺ + C7H2 --> C7H2⁺ + C2H2 ), 
    (@reaction n_H * k219, C2H2⁺ + CH2CCH --> C3H3⁺ + C2H2 ), 
    (@reaction n_H * k220, C2H2⁺ + CH3CCH --> C3H4⁺ + C2H2 ), 
    (@reaction n_H * k221, C2H2⁺ + Fe --> Fe⁺ + C2H2 ), 
    (@reaction n_H * k222, C2H2⁺ + H2CO --> H2CO⁺ + C2H2 ), 
    (@reaction n_H * k223, C2H2⁺ + H2S --> H2S⁺ + C2H2 ), 
    (@reaction n_H * k224, C2H2⁺ + HC4H --> C4H2⁺ + C2H2 ), 
    (@reaction n_H * k225, C2H2⁺ + HCO --> HCO⁺ + C2H2 ), 
    (@reaction n_H * k226, C2H2⁺ + NO --> NO⁺ + C2H2 ), 
    (@reaction n_H * k227, C2H2 + CH2CCH⁺ --> C3H3⁺ + C2H2 ), 
    (@reaction n_H * k228, C2H2 + CO2⁺ --> CO2 + C2H2⁺ ), 
    (@reaction n_H * k229, C2H2 + HC3N⁺ --> C2H2⁺ + HC3N ), 
    (@reaction n_H * k230, C2H2 + HCN⁺ --> C2H2⁺ + HCN ), 
    (@reaction n_H * k231, C2H3 + C2H4⁺ --> C2H4 + C2H3⁺ ), 
    (@reaction n_H * k232, C2H4⁺ + CH2CCH --> C3H3⁺ + C2H4 ), 
    (@reaction n_H * k233, C2H4⁺ + H2S --> H2S⁺ + C2H4 ), 
    (@reaction n_H * k234, C2H4 + C2N2⁺ --> C2H4⁺ + NCCN ), 
    (@reaction n_H * k235, C2H4 + C3⁺ --> C2H4⁺ + C3 ), 
    (@reaction n_H * k236, C2H4 + C5⁺ --> C2H4⁺ + C5 ), 
    (@reaction n_H * k237, C2H4 + CH3CH3⁺ --> CH3CH3 + C2H4⁺ ), 
    (@reaction n_H * k238, C2H4 + CO2⁺ --> CO2 + C2H4⁺ ), 
    (@reaction n_H * k239, C2H4 + HC3N⁺ --> C2H4⁺ + HC3N ), 
    (@reaction n_H * k240, C2H4 + O2⁺ --> C2H4⁺ + O2 ), 
    (@reaction n_H * k241, C2H + CN⁺ --> CN + C2H⁺ ), 
    (@reaction n_H * k242, C2H + CO⁺ --> CO + C2H⁺ ), 
    (@reaction n_H * k243, C2H + N2⁺ --> N2 + C2H⁺ ), 
    (@reaction n_H * k244, C + C2⁺ --> C2 + C⁺ ), 
    (@reaction n_H * k245, C + CN⁺ --> CN + C⁺ ), 
    (@reaction n_H * k246, C + CO⁺ --> CO + C⁺ ), 
    (@reaction n_H * k247, C + N2⁺ --> N2 + C⁺ ), 
    (@reaction n_H * k248, C + O2⁺ --> O2 + C⁺ ), 
    (@reaction n_H * k249, C + PN⁺ --> PN + C⁺ ), 
    (@reaction n_H * k250, CH⁺ + Fe --> Fe⁺ + CH ), 
    (@reaction n_H * k251, CH⁺ + HCO --> HCO⁺ + CH ), 
    (@reaction n_H * k252, CH⁺ + Mg --> Mg⁺ + CH ), 
    (@reaction n_H * k253, CH⁺ + NH3 --> NH3⁺ + CH ), 
    (@reaction n_H * k254, CH⁺ + NO --> NO⁺ + CH ), 
    (@reaction n_H * k255, CH⁺ + Na --> Na⁺ + CH ), 
    (@reaction n_H * k256, CH⁺ + S --> S⁺ + CH ), 
    (@reaction n_H * k257, CH⁺ + Si --> Si⁺ + CH ), 
    (@reaction n_H * k258, CH2⁺ + NO --> NO⁺ + CH2 ), 
    (@reaction n_H * k259, CH2 + C2⁺ --> C2 + CH2⁺ ), 
    (@reaction n_H * k260, CH2 + CN⁺ --> CN + CH2⁺ ), 
    (@reaction n_H * k261, CH2 + CO⁺ --> CO + CH2⁺ ), 
    (@reaction n_H * k262, CH2 + H2CO⁺ --> H2CO + CH2⁺ ), 
    (@reaction n_H * k263, CH2 + H2O⁺ --> H2O + CH2⁺ ), 
    (@reaction n_H * k264, CH2 + N2⁺ --> N2 + CH2⁺ ), 
    (@reaction n_H * k265, CH2 + NH2⁺ --> NH2 + CH2⁺ ), 
    (@reaction n_H * k266, CH2 + O⁺ --> O + CH2⁺ ), 
    (@reaction n_H * k267, CH2 + O2⁺ --> O2 + CH2⁺ ), 
    (@reaction n_H * k268, CH2 + OH⁺ --> OH + CH2⁺ ), 
    (@reaction n_H * k269, CH3⁺ + C2H3 --> C2H3⁺ + CH3 ), 
    (@reaction n_H * k270, CH3⁺ + CH3COCH3 --> CH3COCH3⁺ + CH3 ), 
    (@reaction n_H * k271, CH3⁺ + Fe --> Fe⁺ + CH3 ), 
    (@reaction n_H * k272, CH3⁺ + HCO --> HCO⁺ + CH3 ), 
    (@reaction n_H * k273, CH3⁺ + Mg --> Mg⁺ + CH3 ), 
    (@reaction n_H * k274, CH3⁺ + NO --> NO⁺ + CH3 ), 
    (@reaction n_H * k275, CH3⁺ + Na --> Na⁺ + CH3 ), 
    (@reaction n_H * k276, CH3CH3⁺ + H2S --> H2S⁺ + CH3CH3 ), 
    (@reaction n_H * k277, CH4⁺ + C2H2 --> C2H2⁺ + CH4 ), 
    (@reaction n_H * k278, CH4⁺ + C2H4 --> C2H4⁺ + CH4 ), 
    (@reaction n_H * k279, CH4⁺ + CH3OH --> CH3OH⁺ + CH4 ), 
    (@reaction n_H * k280, CH4⁺ + H2CO --> H2CO⁺ + CH4 ), 
    (@reaction n_H * k281, CH4⁺ + H2S --> H2S⁺ + CH4 ), 
    (@reaction n_H * k282, CH4⁺ + NH3 --> NH3⁺ + CH4 ), 
    (@reaction n_H * k283, CH4⁺ + O2 --> O2⁺ + CH4 ), 
    (@reaction n_H * k284, CH4⁺ + OCS --> OCS⁺ + CH4 ), 
    (@reaction n_H * k285, CH4 + CO⁺ --> CO + CH4⁺ ), 
    (@reaction n_H * k286, CH4 + CO2⁺ --> CO2 + CH4⁺ ), 
    (@reaction n_H * k287, CH + C2⁺ --> C2 + CH⁺ ), 
    (@reaction n_H * k288, CH + CN⁺ --> CN + CH⁺ ), 
    (@reaction n_H * k289, CH + CO⁺ --> CO + CH⁺ ), 
    (@reaction n_H * k290, CH + H2CO⁺ --> H2CO + CH⁺ ), 
    (@reaction n_H * k291, CH + H2O⁺ --> H2O + CH⁺ ), 
    (@reaction n_H * k292, CH + N⁺ --> N + CH⁺ ), 
    (@reaction n_H * k293, CH + N2⁺ --> N2 + CH⁺ ), 
    (@reaction n_H * k294, CH + NH2⁺ --> NH2 + CH⁺ ), 
    (@reaction n_H * k295, CH + O⁺ --> O + CH⁺ ), 
    (@reaction n_H * k296, CH + O2⁺ --> O2 + CH⁺ ), 
    (@reaction n_H * k297, CH + OH⁺ --> OH + CH⁺ ), 
    (@reaction n_H * k298, CN⁺ + CO2 --> CO2⁺ + CN ), 
    (@reaction n_H * k299, CN⁺ + CO --> CO⁺ + CN ), 
    (@reaction n_H * k300, CN⁺ + H2CO --> H2CO⁺ + CN ), 
    (@reaction n_H * k301, CN⁺ + HCN --> HCN⁺ + CN ), 
    (@reaction n_H * k302, CN⁺ + HCO --> HCO⁺ + CN ), 
    (@reaction n_H * k303, CN⁺ + NO --> NO⁺ + CN ), 
    (@reaction n_H * k304, CN⁺ + O2 --> O2⁺ + CN ), 
    (@reaction n_H * k305, CN⁺ + S --> S⁺ + CN ), 
    (@reaction n_H * k306, CN + N2⁺ --> N2 + CN⁺ ), 
    (@reaction n_H * k307, CO⁺ + CO2 --> CO2⁺ + CO ), 
    (@reaction n_H * k308, CO⁺ + H2CO --> H2CO⁺ + CO ), 
    (@reaction n_H * k309, CO⁺ + H2S --> H2S⁺ + CO ), 
    (@reaction n_H * k310, CO⁺ + HCO --> HCO⁺ + CO ), 
    (@reaction n_H * k311, CO⁺ + NO --> NO⁺ + CO ), 
    (@reaction n_H * k312, CO⁺ + O2 --> O2⁺ + CO ), 
    (@reaction n_H * k313, CO⁺ + S --> S⁺ + CO ), 
    (@reaction n_H * k314, CO2⁺ + OCS --> OCS⁺ + CO2 ), 
    (@reaction n_H * k315, CO2⁺ + SO2 --> SO2⁺ + CO2 ), 
    (@reaction n_H * k316, CO + N2⁺ --> N2 + CO⁺ ), 
    (@reaction n_H * k317, CS⁺ + Fe --> Fe⁺ + CS ), 
    (@reaction n_H * k318, Cl⁺ + H --> Cl + H⁺ ), 
    (@reaction n_H * k319, Cl + H⁺ --> Cl⁺ + H ), 
    (@reaction n_H * k320, H⁺ + C10 --> C10⁺ + H ), 
    (@reaction n_H * k321, H⁺ + C10H2 --> C10H2⁺ + H ), 
    (@reaction n_H * k322, H⁺ + C10H --> C10H⁺ + H ), 
    (@reaction n_H * k323, H⁺ + C2 --> C2⁺ + H ), 
    (@reaction n_H * k324, H⁺ + C2H2 --> C2H2⁺ + H ), 
    (@reaction n_H * k325, H⁺ + C2H3 --> C2H3⁺ + H ), 
    (@reaction n_H * k326, H⁺ + C2H4 --> C2H4⁺ + H ), 
    (@reaction n_H * k327, H⁺ + C2H5 --> C2H3⁺ + H2 + H ), 
    (@reaction n_H * k328, H⁺ + C2H5OH --> C2H5OH⁺ + H ), 
    (@reaction n_H * k329, H⁺ + C2H --> C2H⁺ + H ), 
    (@reaction n_H * k330, H⁺ + C2N --> C2N⁺ + H ), 
    (@reaction n_H * k331, H⁺ + C2O --> C2O⁺ + H ), 
    (@reaction n_H * k332, H⁺ + C2S --> C2S⁺ + H ), 
    (@reaction n_H * k333, H⁺ + C3 --> C3⁺ + H ), 
    (@reaction n_H * k334, H⁺ + C3H2 --> C3H2⁺ + H ), 
    (@reaction n_H * k335, H⁺ + C3H --> C3H⁺ + H ), 
    (@reaction n_H * k336, H⁺ + C3O --> C3O⁺ + H ), 
    (@reaction n_H * k337, H⁺ + C3S --> C3S⁺ + H ), 
    (@reaction n_H * k338, H⁺ + C4 --> C4⁺ + H ), 
    (@reaction n_H * k339, H⁺ + C4H3 --> C4H3⁺ + H ), 
    (@reaction n_H * k340, H⁺ + C4H --> C4H⁺ + H ), 
    (@reaction n_H * k341, H⁺ + C4P --> C4P⁺ + H ), 
    (@reaction n_H * k342, H⁺ + C4S --> C4S⁺ + H ), 
    (@reaction n_H * k343, H⁺ + C5 --> C5⁺ + H ), 
    (@reaction n_H * k344, H⁺ + C5H2 --> C5H2⁺ + H ), 
    (@reaction n_H * k345, H⁺ + C5H --> C5H⁺ + H ), 
    (@reaction n_H * k346, H⁺ + C6 --> C6⁺ + H ), 
    (@reaction n_H * k347, H⁺ + C6H2 --> C6H2⁺ + H ), 
    (@reaction n_H * k348, H⁺ + C6H --> C6H⁺ + H ), 
    (@reaction n_H * k349, H⁺ + C7 --> C7⁺ + H ), 
    (@reaction n_H * k350, H⁺ + C7H2 --> C7H2⁺ + H ), 
    (@reaction n_H * k351, H⁺ + C7H --> C7H⁺ + H ), 
    (@reaction n_H * k352, H⁺ + C8 --> C8⁺ + H ), 
    (@reaction n_H * k353, H⁺ + C8H2 --> C8H2⁺ + H ), 
    (@reaction n_H * k354, H⁺ + C8H --> C8H⁺ + H ), 
    (@reaction n_H * k355, H⁺ + C9 --> C9⁺ + H ), 
    (@reaction n_H * k356, H⁺ + C9H2 --> C9H2⁺ + H ), 
    (@reaction n_H * k357, H⁺ + C9H --> C9H⁺ + H ), 
    (@reaction n_H * k358, H⁺ + CCP --> CCP⁺ + H ), 
    (@reaction n_H * k359, H⁺ + CH2 --> CH2⁺ + H ), 
    (@reaction n_H * k360, H⁺ + CH2CCH --> C3H3⁺ + H ), 
    (@reaction n_H * k361, H⁺ + CH2CN --> CH2CN⁺ + H ), 
    (@reaction n_H * k362, H⁺ + CH3 --> CH3⁺ + H ), 
    (@reaction n_H * k363, H⁺ + CH3C4H --> CH3C4H⁺ + H ), 
    (@reaction n_H * k364, H⁺ + CH3C6H --> C7H4⁺ + H ), 
    (@reaction n_H * k365, H⁺ + CH3CCH --> C3H4⁺ + H ), 
    (@reaction n_H * k366, H⁺ + CH3CHO --> CH3CHO⁺ + H ), 
    (@reaction n_H * k367, H⁺ + CH3CN --> CH3CN⁺ + H ), 
    (@reaction n_H * k368, H⁺ + CH3COCH3 --> CH3COCH3⁺ + H ), 
    (@reaction n_H * k369, H⁺ + CH3OCH3 --> CH3OCH3⁺ + H ), 
    (@reaction n_H * k370, H⁺ + CH3OH --> CH3OH⁺ + H ), 
    (@reaction n_H * k371, H⁺ + CH4 --> CH4⁺ + H ), 
    (@reaction n_H * k372, H⁺ + CH --> CH⁺ + H ), 
    (@reaction n_H * k373, H⁺ + CP --> CP⁺ + H ), 
    (@reaction n_H * k374, H⁺ + CS --> CS⁺ + H ), 
    (@reaction n_H * k375, H⁺ + Fe --> Fe⁺ + H ), 
    (@reaction n_H * k376, H⁺ + H2CO --> H2CO⁺ + H ), 
    (@reaction n_H * k377, H⁺ + H2CS --> H2CS⁺ + H ), 
    (@reaction n_H * k378, H⁺ + H2O --> H2O⁺ + H ), 
    (@reaction n_H * k379, H⁺ + H2S2 --> H2S2⁺ + H ), 
    (@reaction n_H * k380, H⁺ + H2S --> H2S⁺ + H ), 
    (@reaction n_H * k381, H⁺ + H2SiO --> H2SiO⁺ + H ), 
    (@reaction n_H * k382, H⁺ + HC2P --> HC2P⁺ + H ), 
    (@reaction n_H * k383, H⁺ + HC3N --> HC3N⁺ + H ), 
    (@reaction n_H * k384, H⁺ + HC4H --> C4H2⁺ + H ), 
    (@reaction n_H * k385, H⁺ + HC5N --> HC5N⁺ + H ), 
    (@reaction n_H * k386, H⁺ + HC7N --> HC7N⁺ + H ), 
    (@reaction n_H * k387, H⁺ + HC9N --> HC9N⁺ + H ), 
    (@reaction n_H * k388, H⁺ + HCN --> HCN⁺ + H ), 
    (@reaction n_H * k389, H⁺ + HCO --> HCO⁺ + H ), 
    (@reaction n_H * k390, H⁺ + HCOOCH3 --> COOCH4⁺ + H ), 
    (@reaction n_H * k391, H⁺ + HCP --> HCP⁺ + H ), 
    (@reaction n_H * k392, H⁺ + HCSi --> HCSi⁺ + H ), 
    (@reaction n_H * k393, H⁺ + HCl --> HCl⁺ + H ), 
    (@reaction n_H * k394, H⁺ + HNSi --> HNSi⁺ + H ), 
    (@reaction n_H * k395, H⁺ + HPO --> HPO⁺ + H ), 
    (@reaction n_H * k396, H⁺ + HS2 --> HS2⁺ + H ), 
    (@reaction n_H * k397, H⁺ + HS --> HS⁺ + H ), 
    (@reaction n_H * k398, H⁺ + Mg --> Mg⁺ + H ), 
    (@reaction n_H * k399, H⁺ + N2O --> N2O⁺ + H ), 
    (@reaction n_H * k400, H⁺ + NH2 --> NH2⁺ + H ), 
    (@reaction n_H * k401, H⁺ + NH3 --> NH3⁺ + H ), 
    (@reaction n_H * k402, H⁺ + NH --> NH⁺ + H ), 
    (@reaction n_H * k403, H⁺ + NO --> NO⁺ + H ), 
    (@reaction n_H * k404, H⁺ + NS --> NS⁺ + H ), 
    (@reaction n_H * k405, H⁺ + O2 --> O2⁺ + H ), 
    (@reaction n_H * k406, H⁺ + O --> O⁺ + H ), 
    (@reaction n_H * k407, H⁺ + OCS --> OCS⁺ + H ), 
    (@reaction n_H * k408, H⁺ + OH --> OH⁺ + H ), 
    (@reaction n_H * k409, H⁺ + P --> P⁺ + H ), 
    (@reaction n_H * k410, H⁺ + PH2 --> PH2⁺ + H ), 
    (@reaction n_H * k411, H⁺ + PH --> PH⁺ + H ), 
    (@reaction n_H * k412, H⁺ + PN --> PN⁺ + H ), 
    (@reaction n_H * k413, H⁺ + PO --> PO⁺ + H ), 
    (@reaction n_H * k414, H⁺ + S2 --> S2⁺ + H ), 
    (@reaction n_H * k415, H⁺ + S --> S⁺ + H ), 
    (@reaction n_H * k416, H⁺ + SO2 --> SO2⁺ + H ), 
    (@reaction n_H * k417, H⁺ + SO --> SO⁺ + H ), 
    (@reaction n_H * k418, H⁺ + Si --> Si⁺ + H ), 
    (@reaction n_H * k419, H⁺ + SiC2 --> SiC2⁺ + H ), 
    (@reaction n_H * k420, H⁺ + SiC2H2 --> SiC2H2⁺ + H ), 
    (@reaction n_H * k421, H⁺ + SiC2H --> SiC2H⁺ + H ), 
    (@reaction n_H * k422, H⁺ + SiC3 --> SiC3⁺ + H ), 
    (@reaction n_H * k423, H⁺ + SiC3H --> SiC3H⁺ + H ), 
    (@reaction n_H * k424, H⁺ + SiC4 --> SiC4⁺ + H ), 
    (@reaction n_H * k425, H⁺ + SiC --> SiC⁺ + H ), 
    (@reaction n_H * k426, H⁺ + SiCH2 --> SiCH2⁺ + H ), 
    (@reaction n_H * k427, H⁺ + SiCH3 --> SiCH3⁺ + H ), 
    (@reaction n_H * k428, H⁺ + SiH2 --> SiH2⁺ + H ), 
    (@reaction n_H * k429, H⁺ + SiH3 --> SiH3⁺ + H ), 
    (@reaction n_H * k430, H⁺ + SiH4 --> SiH4⁺ + H ), 
    (@reaction n_H * k431, H⁺ + SiH --> SiH⁺ + H ), 
    (@reaction n_H * k432, H⁺ + SiN --> SiN⁺ + H ), 
    (@reaction n_H * k433, H⁺ + SiNC --> SiNC⁺ + H ), 
    (@reaction n_H * k434, H⁺ + SiO --> SiO⁺ + H ), 
    (@reaction n_H * k435, H⁺ + SiS --> SiS⁺ + H ), 
    (@reaction n_H * k436, H2⁺ + C2 --> C2⁺ + H2 ), 
    (@reaction n_H * k437, H2⁺ + C2H2 --> C2H2⁺ + H2 ), 
    (@reaction n_H * k438, H2⁺ + C2H4 --> C2H4⁺ + H2 ), 
    (@reaction n_H * k439, H2⁺ + C2H --> C2H⁺ + H2 ), 
    (@reaction n_H * k440, H2⁺ + CH2 --> CH2⁺ + H2 ), 
    (@reaction n_H * k441, H2⁺ + CH3CH3 --> CH3CH3⁺ + H2 ), 
    (@reaction n_H * k442, H2⁺ + CH4 --> CH4⁺ + H2 ), 
    (@reaction n_H * k443, H2⁺ + CH --> CH⁺ + H2 ), 
    (@reaction n_H * k444, H2⁺ + CN --> CN⁺ + H2 ), 
    (@reaction n_H * k445, H2⁺ + CO --> CO⁺ + H2 ), 
    (@reaction n_H * k446, H2⁺ + H2CO --> H2CO⁺ + H2 ), 
    (@reaction n_H * k447, H2⁺ + H2O --> H2O⁺ + H2 ), 
    (@reaction n_H * k448, H2⁺ + H2S --> H2S⁺ + H2 ), 
    (@reaction n_H * k449, H2⁺ + HCN --> HCN⁺ + H2 ), 
    (@reaction n_H * k450, H2⁺ + HCO --> HCO⁺ + H2 ), 
    (@reaction n_H * k451, H2⁺ + NH2 --> NH2⁺ + H2 ), 
    (@reaction n_H * k452, H2⁺ + NH3 --> NH3⁺ + H2 ), 
    (@reaction n_H * k453, H2⁺ + NH --> NH⁺ + H2 ), 
    (@reaction n_H * k454, H2⁺ + NO --> NO⁺ + H2 ), 
    (@reaction n_H * k455, H2⁺ + O2 --> O2⁺ + H2 ), 
    (@reaction n_H * k456, H2⁺ + OH --> OH⁺ + H2 ), 
    (@reaction n_H * k457, H2 + F⁺ --> H2⁺ + F ), 
    (@reaction n_H * k458, H2 + He⁺ --> He + H2⁺ ), 
    (@reaction n_H * k459, H2CO⁺ + Fe --> Fe⁺ + H2CO ), 
    (@reaction n_H * k460, H2CO⁺ + S --> S⁺ + H2CO ), 
    (@reaction n_H * k461, H2CO + O2⁺ --> O2 + H2CO⁺ ), 
    (@reaction n_H * k462, H2O⁺ + C2 --> C2⁺ + H2O ), 
    (@reaction n_H * k463, H2O⁺ + C2H2 --> C2H2⁺ + H2O ), 
    (@reaction n_H * k464, H2O⁺ + C2H4 --> C2H4⁺ + H2O ), 
    (@reaction n_H * k465, H2O⁺ + C2H --> C2H⁺ + H2O ), 
    (@reaction n_H * k466, H2O⁺ + CH3CH3 --> CH3CH3⁺ + H2O ), 
    (@reaction n_H * k467, H2O⁺ + Fe --> Fe⁺ + H2O ), 
    (@reaction n_H * k468, H2O⁺ + H2CO --> H2CO⁺ + H2O ), 
    (@reaction n_H * k469, H2O⁺ + H2S --> H2S⁺ + H2O ), 
    (@reaction n_H * k470, H2O⁺ + HCO --> HCO⁺ + H2O ), 
    (@reaction n_H * k471, H2O⁺ + Mg --> Mg⁺ + H2O ), 
    (@reaction n_H * k472, H2O⁺ + NO --> NO⁺ + H2O ), 
    (@reaction n_H * k473, H2O⁺ + Na --> Na⁺ + H2O ), 
    (@reaction n_H * k474, H2O⁺ + O2 --> O2⁺ + H2O ), 
    (@reaction n_H * k475, H2O⁺ + OCS --> OCS⁺ + H2O ), 
    (@reaction n_H * k476, H2O⁺ + S --> S⁺ + H2O ), 
    (@reaction n_H * k477, H2O⁺ + Si --> Si⁺ + H2O ), 
    (@reaction n_H * k478, H2O + CO⁺ --> CO + H2O⁺ ), 
    (@reaction n_H * k479, H2O + CO2⁺ --> CO2 + H2O⁺ ), 
    (@reaction n_H * k480, H2O + HCN⁺ --> HCN + H2O⁺ ), 
    (@reaction n_H * k481, H2O + N2⁺ --> N2 + H2O⁺ ), 
    (@reaction n_H * k482, H2O + N2O⁺ --> H2O⁺ + N2O ), 
    (@reaction n_H * k483, H2S⁺ + Fe --> Fe⁺ + H2S ), 
    (@reaction n_H * k484, H2S + CO2⁺ --> CO2 + H2S⁺ ), 
    (@reaction n_H * k485, H2S + N2O⁺ --> H2S⁺ + N2O ), 
    (@reaction n_H * k486, H2S + OCS⁺ --> H2S⁺ + OCS ), 
    (@reaction n_H * k487, H + CN⁺ --> CN + H⁺ ), 
    (@reaction n_H * k488, H + CO⁺ --> CO + H⁺ ), 
    (@reaction n_H * k489, H + H2⁺ --> H2 + H⁺ ), 
    (@reaction n_H * k490, H + HCN⁺ --> HCN + H⁺ ), 
    (@reaction n_H * k491, H + He⁺ --> He + H⁺ ), 
    (@reaction n_H * k492, H + O⁺ --> O + H⁺ ), 
    (@reaction n_H * k493, HC4H + HC3N⁺ --> C4H2⁺ + HC3N ), 
    (@reaction n_H * k494, HCN⁺ + NO --> NO⁺ + HCN ), 
    (@reaction n_H * k495, HCN⁺ + O2 --> O2⁺ + HCN ), 
    (@reaction n_H * k496, HCN⁺ + S --> S⁺ + HCN ), 
    (@reaction n_H * k497, HCN + CO⁺ --> CO + HCN⁺ ), 
    (@reaction n_H * k498, HCN + N2⁺ --> N2 + HCN⁺ ), 
    (@reaction n_H * k499, HCO⁺ + Fe --> Fe⁺ + HCO ), 
    (@reaction n_H * k500, HCO + H2CO⁺ --> H2CO + HCO⁺ ), 
    (@reaction n_H * k501, HCO + H2S⁺ --> H2S + HCO⁺ ), 
    (@reaction n_H * k502, HCO + O2⁺ --> O2 + HCO⁺ ), 
    (@reaction n_H * k503, HCO + S⁺ --> S + HCO⁺ ), 
    (@reaction n_H * k504, HCO + SiO⁺ --> SiO + HCO⁺ ), 
    (@reaction n_H * k505, HS⁺ + Fe --> Fe⁺ + HS ), 
    (@reaction n_H * k506, He⁺ + C10 --> C10⁺ + He ), 
    (@reaction n_H * k507, He⁺ + C10H --> C10H⁺ + He ), 
    (@reaction n_H * k508, He⁺ + C2 --> C2⁺ + He ), 
    (@reaction n_H * k509, He⁺ + C2H2 --> C2H2⁺ + He ), 
    (@reaction n_H * k510, He⁺ + C2H4 --> C2H4⁺ + He ), 
    (@reaction n_H * k511, He⁺ + C2H5 --> C2H5⁺ + He ), 
    (@reaction n_H * k512, He⁺ + C3N --> C3N⁺ + He ), 
    (@reaction n_H * k513, He⁺ + C --> C⁺ + He ), 
    (@reaction n_H * k514, He⁺ + CH4 --> CH4⁺ + He ), 
    (@reaction n_H * k515, He⁺ + CH --> CH⁺ + He ), 
    (@reaction n_H * k516, He⁺ + CO2 --> CO2⁺ + He ), 
    (@reaction n_H * k517, He⁺ + H2CO --> H2CO⁺ + He ), 
    (@reaction n_H * k518, He⁺ + H2O --> H2O⁺ + He ), 
    (@reaction n_H * k519, He⁺ + H2S --> H2S⁺ + He ), 
    (@reaction n_H * k520, He⁺ + HCNO --> HCNO⁺ + He ), 
    (@reaction n_H * k521, He⁺ + HNCO --> HNCO⁺ + He ), 
    (@reaction n_H * k522, He⁺ + HONC --> HONC⁺ + He ), 
    (@reaction n_H * k523, He⁺ + N2 --> N2⁺ + He ), 
    (@reaction n_H * k524, He⁺ + NH3 --> NH3⁺ + He ), 
    (@reaction n_H * k525, He⁺ + O2 --> O2⁺ + He ), 
    (@reaction n_H * k526, He⁺ + P --> P⁺ + He ), 
    (@reaction n_H * k527, He⁺ + SO2 --> SO2⁺ + He ), 
    (@reaction n_H * k528, He⁺ + Si --> Si⁺ + He ), 
    (@reaction n_H * k529, Mg + C2H2⁺ --> C2H2 + Mg⁺ ), 
    (@reaction n_H * k530, Mg + CS⁺ --> CS + Mg⁺ ), 
    (@reaction n_H * k531, Mg + H2CO⁺ --> H2CO + Mg⁺ ), 
    (@reaction n_H * k532, Mg + H2S⁺ --> H2S + Mg⁺ ), 
    (@reaction n_H * k533, Mg + HCO⁺ --> HCO + Mg⁺ ), 
    (@reaction n_H * k534, Mg + HS⁺ --> HS + Mg⁺ ), 
    (@reaction n_H * k535, Mg + N2⁺ --> N2 + Mg⁺ ), 
    (@reaction n_H * k536, Mg + NO⁺ --> NO + Mg⁺ ), 
    (@reaction n_H * k537, Mg + O2⁺ --> O2 + Mg⁺ ), 
    (@reaction n_H * k538, Mg + S⁺ --> S + Mg⁺ ), 
    (@reaction n_H * k539, Mg + SO⁺ --> SO + Mg⁺ ), 
    (@reaction n_H * k540, Mg + Si⁺ --> Si + Mg⁺ ), 
    (@reaction n_H * k541, Mg + SiO⁺ --> SiO + Mg⁺ ), 
    (@reaction n_H * k542, N⁺ + C2 --> C2⁺ + N ), 
    (@reaction n_H * k543, N⁺ + C2H --> C2H⁺ + N ), 
    (@reaction n_H * k544, N⁺ + CH2 --> CH2⁺ + N ), 
    (@reaction n_H * k545, N⁺ + CH3OH --> CH3OH⁺ + N ), 
    (@reaction n_H * k546, N⁺ + CH4 --> CH4⁺ + N ), 
    (@reaction n_H * k547, N⁺ + CN --> CN⁺ + N ), 
    (@reaction n_H * k548, N⁺ + CO2 --> CO2⁺ + N ), 
    (@reaction n_H * k549, N⁺ + CO --> CO⁺ + N ), 
    (@reaction n_H * k550, N⁺ + Fe --> Fe⁺ + N ), 
    (@reaction n_H * k551, N⁺ + H2CO --> H2CO⁺ + N ), 
    (@reaction n_H * k552, N⁺ + H2O --> H2O⁺ + N ), 
    (@reaction n_H * k553, N⁺ + H2S --> H2S⁺ + N ), 
    (@reaction n_H * k554, N⁺ + HCN --> HCN⁺ + N ), 
    (@reaction n_H * k555, N⁺ + HCO --> HCO⁺ + N ), 
    (@reaction n_H * k556, N⁺ + Mg --> Mg⁺ + N ), 
    (@reaction n_H * k557, N⁺ + NH2 --> NH2⁺ + N ), 
    (@reaction n_H * k558, N⁺ + NH3 --> NH3⁺ + N ), 
    (@reaction n_H * k559, N⁺ + NH --> NH⁺ + N ), 
    (@reaction n_H * k560, N⁺ + NO --> NO⁺ + N ), 
    (@reaction n_H * k561, N⁺ + O2 --> O2⁺ + N ), 
    (@reaction n_H * k562, N⁺ + OCS --> OCS⁺ + N ), 
    (@reaction n_H * k563, N⁺ + OH --> OH⁺ + N ), 
    (@reaction n_H * k564, N2⁺ + CO2 --> CO2⁺ + N2 ), 
    (@reaction n_H * k565, N2⁺ + Fe --> Fe⁺ + N2 ), 
    (@reaction n_H * k566, N2⁺ + H2CO --> H2CO⁺ + N2 ), 
    (@reaction n_H * k567, N2⁺ + H2S --> H2S⁺ + N2 ), 
    (@reaction n_H * k568, N2⁺ + HCO --> HCO⁺ + N2 ), 
    (@reaction n_H * k569, N2⁺ + NO --> NO⁺ + N2 ), 
    (@reaction n_H * k570, N2⁺ + O2 --> O2⁺ + N2 ), 
    (@reaction n_H * k571, N2⁺ + OCS --> OCS⁺ + N2 ), 
    (@reaction n_H * k572, N2⁺ + S --> S⁺ + N2 ), 
    (@reaction n_H * k573, N + N2⁺ --> N2 + N⁺ ), 
    (@reaction n_H * k574, NH⁺ + H2CO --> H2CO⁺ + NH ), 
    (@reaction n_H * k575, NH⁺ + H2O --> H2O⁺ + NH ), 
    (@reaction n_H * k576, NH⁺ + NH3 --> NH3⁺ + NH ), 
    (@reaction n_H * k577, NH⁺ + NO --> NO⁺ + NH ), 
    (@reaction n_H * k578, NH⁺ + O2 --> O2⁺ + NH ), 
    (@reaction n_H * k579, NH⁺ + S --> S⁺ + NH ), 
    (@reaction n_H * k580, NH2⁺ + H2S --> H2S⁺ + NH2 ), 
    (@reaction n_H * k581, NH2⁺ + HCO --> HCO⁺ + NH2 ), 
    (@reaction n_H * k582, NH2⁺ + NH3 --> NH3⁺ + NH2 ), 
    (@reaction n_H * k583, NH2⁺ + NO --> NO⁺ + NH2 ), 
    (@reaction n_H * k584, NH2⁺ + S --> S⁺ + NH2 ), 
    (@reaction n_H * k585, NH2 + C2⁺ --> C2 + NH2⁺ ), 
    (@reaction n_H * k586, NH2 + CN⁺ --> CN + NH2⁺ ), 
    (@reaction n_H * k587, NH2 + CO⁺ --> CO + NH2⁺ ), 
    (@reaction n_H * k588, NH2 + H2O⁺ --> H2O + NH2⁺ ), 
    (@reaction n_H * k589, NH2 + N2⁺ --> N2 + NH2⁺ ), 
    (@reaction n_H * k590, NH2 + O2⁺ --> O2 + NH2⁺ ), 
    (@reaction n_H * k591, NH2 + OH⁺ --> OH + NH2⁺ ), 
    (@reaction n_H * k592, NH3⁺ + Fe --> Fe⁺ + NH3 ), 
    (@reaction n_H * k593, NH3⁺ + HCO --> HCO⁺ + NH3 ), 
    (@reaction n_H * k594, NH3⁺ + Mg --> Mg⁺ + NH3 ), 
    (@reaction n_H * k595, NH3⁺ + NO --> NO⁺ + NH3 ), 
    (@reaction n_H * k596, NH3⁺ + Na --> Na⁺ + NH3 ), 
    (@reaction n_H * k597, NH3⁺ + Si --> Si⁺ + NH3 ), 
    (@reaction n_H * k598, NH3 + C2H2⁺ --> C2H2 + NH3⁺ ), 
    (@reaction n_H * k599, NH3 + C2H4⁺ --> C2H4 + NH3⁺ ), 
    (@reaction n_H * k600, NH3 + C3H⁺ --> C3H + NH3⁺ ), 
    (@reaction n_H * k601, NH3 + CH3CH3⁺ --> CH3CH3 + NH3⁺ ), 
    (@reaction n_H * k602, NH3 + CO⁺ --> CO + NH3⁺ ), 
    (@reaction n_H * k603, NH3 + CO2⁺ --> CO2 + NH3⁺ ), 
    (@reaction n_H * k604, NH3 + H2CO⁺ --> H2CO + NH3⁺ ), 
    (@reaction n_H * k605, NH3 + H2O⁺ --> H2O + NH3⁺ ), 
    (@reaction n_H * k606, NH3 + H2S⁺ --> H2S + NH3⁺ ), 
    (@reaction n_H * k607, NH3 + HC3N⁺ --> HC3N + NH3⁺ ), 
    (@reaction n_H * k608, NH3 + HCN⁺ --> HCN + NH3⁺ ), 
    (@reaction n_H * k609, NH3 + HS⁺ --> HS + NH3⁺ ), 
    (@reaction n_H * k610, NH3 + N2⁺ --> N2 + NH3⁺ ), 
    (@reaction n_H * k611, NH3 + O2⁺ --> O2 + NH3⁺ ), 
    (@reaction n_H * k612, NH3 + OCS⁺ --> NH3⁺ + OCS ), 
    (@reaction n_H * k613, NH3 + P⁺ --> P + NH3⁺ ), 
    (@reaction n_H * k614, NH3 + S⁺ --> S + NH3⁺ ), 
    (@reaction n_H * k615, NH3 + SO⁺ --> SO + NH3⁺ ), 
    (@reaction n_H * k616, NH + CN⁺ --> CN + NH⁺ ), 
    (@reaction n_H * k617, NH + CO⁺ --> CO + NH⁺ ), 
    (@reaction n_H * k618, NH + N2⁺ --> N2 + NH⁺ ), 
    (@reaction n_H * k619, NH + O⁺ --> O + NH⁺ ), 
    (@reaction n_H * k620, NO⁺ + Fe --> Fe⁺ + NO ), 
    (@reaction n_H * k621, NO + C3H⁺ --> C3H + NO⁺ ), 
    (@reaction n_H * k622, NO + C4H2⁺ --> NO⁺ + HC4H ), 
    (@reaction n_H * k623, NO + CO2⁺ --> CO2 + NO⁺ ), 
    (@reaction n_H * k624, NO + H2CO⁺ --> H2CO + NO⁺ ), 
    (@reaction n_H * k625, NO + H2S⁺ --> H2S + NO⁺ ), 
    (@reaction n_H * k626, NO + HNO⁺ --> HNO + NO⁺ ), 
    (@reaction n_H * k627, NO + HS⁺ --> HS + NO⁺ ), 
    (@reaction n_H * k628, NO + O2⁺ --> O2 + NO⁺ ), 
    (@reaction n_H * k629, NO + S⁺ --> S + NO⁺ ), 
    (@reaction n_H * k630, NO + S2⁺ --> S2 + NO⁺ ), 
    (@reaction n_H * k631, NO + SiO⁺ --> SiO + NO⁺ ), 
    (@reaction n_H * k632, Na + C2H2⁺ --> C2H2 + Na⁺ ), 
    (@reaction n_H * k633, Na + CS⁺ --> CS + Na⁺ ), 
    (@reaction n_H * k634, Na + Fe⁺ --> Fe + Na⁺ ), 
    (@reaction n_H * k635, Na + H2CO⁺ --> H2CO + Na⁺ ), 
    (@reaction n_H * k636, Na + H2S⁺ --> H2S + Na⁺ ), 
    (@reaction n_H * k637, Na + HCO⁺ --> HCO + Na⁺ ), 
    (@reaction n_H * k638, Na + HS⁺ --> HS + Na⁺ ), 
    (@reaction n_H * k639, Na + Mg⁺ --> Mg + Na⁺ ), 
    (@reaction n_H * k640, Na + N2⁺ --> N2 + Na⁺ ), 
    (@reaction n_H * k641, Na + NO⁺ --> NO + Na⁺ ), 
    (@reaction n_H * k642, Na + O2⁺ --> O2 + Na⁺ ), 
    (@reaction n_H * k643, Na + S⁺ --> S + Na⁺ ), 
    (@reaction n_H * k644, Na + SO⁺ --> SO + Na⁺ ), 
    (@reaction n_H * k645, Na + Si⁺ --> Si + Na⁺ ), 
    (@reaction n_H * k646, O⁺ + C2 --> C2⁺ + O ), 
    (@reaction n_H * k647, O⁺ + C2H2 --> C2H2⁺ + O ), 
    (@reaction n_H * k648, O⁺ + C2H4 --> C2H4⁺ + O ), 
    (@reaction n_H * k649, O⁺ + C2H --> C2H⁺ + O ), 
    (@reaction n_H * k650, O⁺ + CH3CN --> CH3CN⁺ + O ), 
    (@reaction n_H * k651, O⁺ + CH3OH --> CH3OH⁺ + O ), 
    (@reaction n_H * k652, O⁺ + CH4 --> CH4⁺ + O ), 
    (@reaction n_H * k653, O⁺ + CO --> CO⁺ + O ), 
    (@reaction n_H * k654, O⁺ + Fe --> Fe⁺ + O ), 
    (@reaction n_H * k655, O⁺ + H2CO --> H2CO⁺ + O ), 
    (@reaction n_H * k656, O⁺ + H2O --> H2O⁺ + O ), 
    (@reaction n_H * k657, O⁺ + H2S --> H2S⁺ + O ), 
    (@reaction n_H * k658, O⁺ + HCO --> HCO⁺ + O ), 
    (@reaction n_H * k659, O⁺ + N2O --> N2O⁺ + O ), 
    (@reaction n_H * k660, O⁺ + NH2 --> NH2⁺ + O ), 
    (@reaction n_H * k661, O⁺ + NH3 --> NH3⁺ + O ), 
    (@reaction n_H * k662, O⁺ + O2 --> O2⁺ + O ), 
    (@reaction n_H * k663, O⁺ + OCS --> OCS⁺ + O ), 
    (@reaction n_H * k664, O⁺ + OH --> OH⁺ + O ), 
    (@reaction n_H * k665, O⁺ + SO2 --> SO2⁺ + O ), 
    (@reaction n_H * k666, O⁻ + CN --> CN⁻ + O ), 
    (@reaction n_H * k667, O⁻ + O2 --> O2⁻ + O ), 
    (@reaction n_H * k668, O2⁺ + C2H2 --> C2H2⁺ + O2 ), 
    (@reaction n_H * k669, O2⁺ + CH2CCH2 --> C3H4⁺ + O2 ), 
    (@reaction n_H * k670, O2⁺ + CH3CHCH2 --> C3H6⁺ + O2 ), 
    (@reaction n_H * k671, O2⁺ + CH3OCH3 --> CH3OCH3⁺ + O2 ), 
    (@reaction n_H * k672, O2⁺ + CH3OH --> CH3OH⁺ + O2 ), 
    (@reaction n_H * k673, O2⁺ + Fe --> Fe⁺ + O2 ), 
    (@reaction n_H * k674, O2⁺ + H2S --> H2S⁺ + O2 ), 
    (@reaction n_H * k675, O2⁺ + HCOOH --> HCOOH⁺ + O2 ), 
    (@reaction n_H * k676, O2⁺ + NO2 --> NO2⁺ + O2 ), 
    (@reaction n_H * k677, O2⁺ + S --> S⁺ + O2 ), 
    (@reaction n_H * k678, O2 + CO2⁺ --> CO2 + O2⁺ ), 
    (@reaction n_H * k679, O2 + Cl⁺ --> Cl + O2⁺ ), 
    (@reaction n_H * k680, O2 + SO2⁺ --> SO2 + O2⁺ ), 
    (@reaction n_H * k681, O + CN⁺ --> CN + O⁺ ), 
    (@reaction n_H * k682, O + CO⁺ --> CO + O⁺ ), 
    (@reaction n_H * k683, O + CO2⁺ --> CO2 + O⁺ ), 
    (@reaction n_H * k684, O + N2⁺ --> N2 + O⁺ ), 
    (@reaction n_H * k685, OH⁺ + C2 --> C2⁺ + OH ), 
    (@reaction n_H * k686, OH⁺ + C2H --> C2H⁺ + OH ), 
    (@reaction n_H * k687, OH⁺ + CH3CH3 --> CH3CH3⁺ + OH ), 
    (@reaction n_H * k688, OH⁺ + H2CO --> H2CO⁺ + OH ), 
    (@reaction n_H * k689, OH⁺ + H2O --> H2O⁺ + OH ), 
    (@reaction n_H * k690, OH⁺ + H2S --> H2S⁺ + OH ), 
    (@reaction n_H * k691, OH⁺ + HCO --> HCO⁺ + OH ), 
    (@reaction n_H * k692, OH⁺ + NH3 --> NH3⁺ + OH ), 
    (@reaction n_H * k693, OH⁺ + NO --> NO⁺ + OH ), 
    (@reaction n_H * k694, OH⁺ + O2 --> O2⁺ + OH ), 
    (@reaction n_H * k695, OH⁺ + S --> S⁺ + OH ), 
    (@reaction n_H * k696, OH⁻ + CN --> CN⁻ + OH ), 
    (@reaction n_H * k697, OH + C2⁺ --> C2 + OH⁺ ), 
    (@reaction n_H * k698, OH + CN⁺ --> CN + OH⁺ ), 
    (@reaction n_H * k699, OH + CO⁺ --> CO + OH⁺ ), 
    (@reaction n_H * k700, OH + N2⁺ --> N2 + OH⁺ ), 
    (@reaction n_H * k701, P⁺ + H2S --> H2S⁺ + P ), 
    (@reaction n_H * k702, S⁺ + CH3CHCH2 --> C3H6⁺ + S ), 
    (@reaction n_H * k703, S⁺ + Fe --> Fe⁺ + S ), 
    (@reaction n_H * k704, S⁺ + HC4H --> C4H2⁺ + S ), 
    (@reaction n_H * k705, S⁺ + SiC --> SiC⁺ + S ), 
    (@reaction n_H * k706, S⁺ + SiS --> SiS⁺ + S ), 
    (@reaction n_H * k707, S + C⁺ --> C + S⁺ ), 
    (@reaction n_H * k708, S + H2S⁺ --> H2S + S⁺ ), 
    (@reaction n_H * k709, S + HS⁺ --> HS + S⁺ ), 
    (@reaction n_H * k710, SO⁺ + CH2CCH2 --> C3H4⁺ + SO ), 
    (@reaction n_H * k711, SO⁺ + CH3CHCH2 --> C3H6⁺ + SO ), 
    (@reaction n_H * k712, SO⁺ + CH3COCH3 --> CH3COCH3⁺ + SO ), 
    (@reaction n_H * k713, SO⁺ + CH3OCH3 --> CH3OCH3⁺ + SO ), 
    (@reaction n_H * k714, SO⁺ + Fe --> Fe⁺ + SO ), 
    (@reaction n_H * k715, Si⁺ + Fe --> Fe⁺ + Si ), 
    (@reaction n_H * k716, Si + CS⁺ --> CS + Si⁺ ), 
    (@reaction n_H * k717, Si + H2CO⁺ --> H2CO + Si⁺ ), 
    (@reaction n_H * k718, Si + H2S⁺ --> H2S + Si⁺ ), 
    (@reaction n_H * k719, Si + HS⁺ --> HS + Si⁺ ), 
    (@reaction n_H * k720, Si + NO⁺ --> NO + Si⁺ ), 
    (@reaction n_H * k721, Si + O2⁺ --> O2 + Si⁺ ), 
    (@reaction n_H * k722, Si + P⁺ --> P + Si⁺ ), 
    (@reaction n_H * k723, Si + S⁺ --> S + Si⁺ ), 
    (@reaction n_H * k724, SiH + S⁺ --> S + SiH⁺ ), 
    (@reaction n_H * k725, SiO⁺ + Fe --> Fe⁺ + SiO ), 
    (@reaction n_H * k726, C --> C⁺ + e ), 
    (@reaction n_H * k727, CO --> CO⁺ + e ), 
    (@reaction n_H * k728, Cl --> Cl⁺ + e ), 
    (@reaction n_H * k729, H2 --> H⁺ + H⁻ ), 
    (@reaction n_H * k730, H2 --> H⁺ + H + e ), 
    (@reaction n_H * k731, H2 --> H2⁺ + e ), 
    (@reaction n_H * k732, H2 --> H + H ), 
    (@reaction n_H * k733, H --> H⁺ + e ), 
    (@reaction n_H * k734, He --> He⁺ + e ), 
    (@reaction n_H * k735, N --> N⁺ + e ), 
    (@reaction n_H * k736, O --> O⁺ + e ), 
    (@reaction n_H * k737, C⁻ --> C + e ), 
    (@reaction n_H * k738, C10⁻ --> C10 + e ), 
    (@reaction n_H * k739, C10 --> C9 + C ), 
    (@reaction n_H * k740, C10H⁻ --> C10H + e ), 
    (@reaction n_H * k741, C10H2 --> C10H + H ), 
    (@reaction n_H * k742, C10H --> C10 + H ), 
    (@reaction n_H * k743, C11 --> C10 + C ), 
    (@reaction n_H * k744, C2⁻ --> C2 + e ), 
    (@reaction n_H * k745, C2 --> C + C ), 
    (@reaction n_H * k746, C2H⁻ --> C2H + e ), 
    (@reaction n_H * k747, C2H2 --> C2H2⁺ + e ), 
    (@reaction n_H * k748, C2H2 --> C2H + H ), 
    (@reaction n_H * k749, C2H3 --> C2H2 + H ), 
    (@reaction n_H * k750, C2H3 --> C2H3⁺ + e ), 
    (@reaction n_H * k751, C2H4 --> C2H2 + H2 ), 
    (@reaction n_H * k752, C2H4 --> C2H4⁺ + e ), 
    (@reaction n_H * k753, C2H5 --> C2H3 + H2 ), 
    (@reaction n_H * k754, C2H5 --> C2H5⁺ + e ), 
    (@reaction n_H * k755, C2H5CN --> CH2CHCNH⁺ + H + e ), 
    (@reaction n_H * k756, C2H5CN --> CN + C2H5 ), 
    (@reaction n_H * k757, C2H5OH --> C2H5 + OH ), 
    (@reaction n_H * k758, C2H5OH --> C2H5OH⁺ + e ), 
    (@reaction n_H * k759, C2H --> C2 + H ), 
    (@reaction n_H * k760, C2H --> C2H⁺ + e ), 
    (@reaction n_H * k761, C2N --> C2 + N ), 
    (@reaction n_H * k762, C2N --> CN + C ), 
    (@reaction n_H * k763, C2O --> C2 + O ), 
    (@reaction n_H * k764, C2O --> CO + C ), 
    (@reaction n_H * k765, C2S --> CS + C ), 
    (@reaction n_H * k766, C2S --> S + C2 ), 
    (@reaction n_H * k767, C3⁻ --> C3 + e ), 
    (@reaction n_H * k768, C3 --> C2 + C ), 
    (@reaction n_H * k769, C3H⁻ --> C3H + e ), 
    (@reaction n_H * k770, C3H2 --> C3H + H ), 
    (@reaction n_H * k771, C3H --> C3 + H ), 
    (@reaction n_H * k772, C3N⁻ --> C3N + e ), 
    (@reaction n_H * k773, C3N --> CN + C2 ), 
    (@reaction n_H * k774, C3O --> CO + C2 ), 
    (@reaction n_H * k775, C3P --> CCP + C ), 
    (@reaction n_H * k776, C3S --> CS + C2 ), 
    (@reaction n_H * k777, C4⁻ --> C4 + e ), 
    (@reaction n_H * k778, C4 --> C2 + C2 ), 
    (@reaction n_H * k779, C4 --> C3 + C ), 
    (@reaction n_H * k780, C4H⁻ --> C4H + e ), 
    (@reaction n_H * k781, C4H3 --> C4H + H + H ), 
    (@reaction n_H * k782, C4H --> C2H + C2 ), 
    (@reaction n_H * k783, C4H --> C4 + H ), 
    (@reaction n_H * k784, C4N --> C3 + CN ), 
    (@reaction n_H * k785, C4P --> C3P + C ), 
    (@reaction n_H * k786, C4S --> CS + C3 ), 
    (@reaction n_H * k787, C5⁻ --> C5 + e ), 
    (@reaction n_H * k788, C5 --> C3 + C2 ), 
    (@reaction n_H * k789, C5 --> C4 + C ), 
    (@reaction n_H * k790, C5H⁻ --> C5H + e ), 
    (@reaction n_H * k791, C5H2 --> C3H + C2H ), 
    (@reaction n_H * k792, C5H2 --> C5H + H ), 
    (@reaction n_H * k793, C5H --> C3 + C2H ), 
    (@reaction n_H * k794, C5H --> C3H + C2 ), 
    (@reaction n_H * k795, C5H --> C5 + H ), 
    (@reaction n_H * k796, C5N⁻ --> C5N + e ), 
    (@reaction n_H * k797, C5N --> C4 + CN ), 
    (@reaction n_H * k798, C6⁻ --> C6 + e ), 
    (@reaction n_H * k799, C6 --> C5 + C ), 
    (@reaction n_H * k800, C6H⁻ --> C6H + e ), 
    (@reaction n_H * k801, C6H2 --> C4H + C2H ), 
    (@reaction n_H * k802, C6H6 --> C6H2 + H2 + H + H ), 
    (@reaction n_H * k803, C6H --> C3H + C3 ), 
    (@reaction n_H * k804, C6H --> C4 + C2H ), 
    (@reaction n_H * k805, C6H --> C6 + H ), 
    (@reaction n_H * k806, C7⁻ --> C7 + e ), 
    (@reaction n_H * k807, C7 --> C6 + C ), 
    (@reaction n_H * k808, C7H⁻ --> C7H + e ), 
    (@reaction n_H * k809, C7H2 --> C7H + H ), 
    (@reaction n_H * k810, C7H --> C7 + H ), 
    (@reaction n_H * k811, C7N --> C6 + CN ), 
    (@reaction n_H * k812, C8⁻ --> C8 + e ), 
    (@reaction n_H * k813, C8 --> C7 + C ), 
    (@reaction n_H * k814, C8H⁻ --> C8H + e ), 
    (@reaction n_H * k815, C8H2 --> C8H + H ), 
    (@reaction n_H * k816, C8H --> C8 + H ), 
    (@reaction n_H * k817, C9⁻ --> C9 + e ), 
    (@reaction n_H * k818, C9 --> C8 + C ), 
    (@reaction n_H * k819, C9H⁻ --> C9H + e ), 
    (@reaction n_H * k820, C9H2 --> C9H + H ), 
    (@reaction n_H * k821, C9H --> C9 + H ), 
    (@reaction n_H * k822, C9N --> C8 + CN ), 
    (@reaction n_H * k823, C --> C⁺ + e ), 
    (@reaction n_H * k824, CCP --> C2 + P ), 
    (@reaction n_H * k825, CCP --> CP + C ), 
    (@reaction n_H * k826, CCl --> Cl + C ), 
    (@reaction n_H * k827, CH⁺ --> C⁺ + H ), 
    (@reaction n_H * k828, CH⁻ --> CH + e ), 
    (@reaction n_H * k829, CH2 --> CH2⁺ + e ), 
    (@reaction n_H * k830, CH2 --> CH + H ), 
    (@reaction n_H * k831, CH2CCH2 --> C3H4⁺ + e ), 
    (@reaction n_H * k832, CH2CCH2 --> CH2CCH + H ), 
    (@reaction n_H * k833, CH2CCH --> C3H2 + H ), 
    (@reaction n_H * k834, CH2CCH --> C3H + H2 ), 
    (@reaction n_H * k835, CH2CHCCH --> CH2CCH2 + C ), 
    (@reaction n_H * k836, CH2CHCCH --> CH3CCH + C ), 
    (@reaction n_H * k837, CH2CHCHCH2 --> CH3CHCH2 + C ), 
    (@reaction n_H * k838, CH2CHCN --> C2H3 + CN ), 
    (@reaction n_H * k839, CH2CN --> CH2 + CN ), 
    (@reaction n_H * k840, CH2CO --> CH2CO⁺ + e ), 
    (@reaction n_H * k841, CH2CO --> CO + CH2 ), 
    (@reaction n_H * k842, CH2NH --> NH + CH2 ), 
    (@reaction n_H * k843, CH2PH --> HCP + H2 ), 
    (@reaction n_H * k844, CH3 --> CH2 + H ), 
    (@reaction n_H * k845, CH3 --> CH3⁺ + e ), 
    (@reaction n_H * k846, CH3 --> CH + H2 ), 
    (@reaction n_H * k847, CH3C3N --> C3N + CH3 ), 
    (@reaction n_H * k848, CH3C4H --> C4H + CH3 ), 
    (@reaction n_H * k849, CH3C5N --> CH3 + C5N ), 
    (@reaction n_H * k850, CH3C6H --> CH3 + C6H ), 
    (@reaction n_H * k851, CH3C7N --> CH3 + C7N ), 
    (@reaction n_H * k852, CH3CCH --> C3H4⁺ + e ), 
    (@reaction n_H * k853, CH3CCH --> CH2CCH + H ), 
    (@reaction n_H * k854, CH3CH3 --> C2H4 + H2 ), 
    (@reaction n_H * k855, CH3CH3 --> CH3CH3⁺ + e ), 
    (@reaction n_H * k856, CH3CHCH2 --> C2H4 + CH2 ), 
    (@reaction n_H * k857, CH3CHO --> CH3CHO⁺ + e ), 
    (@reaction n_H * k858, CH3CHO --> CO + CH4 ), 
    (@reaction n_H * k859, CH3CHO --> HCO + CH3 ), 
    (@reaction n_H * k860, CH3CN --> CH3CN⁺ + e ), 
    (@reaction n_H * k861, CH3CN --> CN + CH3 ), 
    (@reaction n_H * k862, CH3COCH3 --> CH2CO + CH4 ), 
    (@reaction n_H * k863, CH3COCH3 --> CH3COCH3⁺ + e ), 
    (@reaction n_H * k864, CH3COCH3 --> CO + CH3 + CH3 ), 
    (@reaction n_H * k865, CH3OCH3 --> CH3OCH3⁺ + e ), 
    (@reaction n_H * k866, CH3OCH3 --> H2CO + CH4 ), 
    (@reaction n_H * k867, CH3OH --> CH3OH⁺ + e ), 
    (@reaction n_H * k868, CH3OH --> H2CO + H2 ), 
    (@reaction n_H * k869, CH3OH --> OH + CH3 ), 
    (@reaction n_H * k870, CH4 --> CH2 + H2 ), 
    (@reaction n_H * k871, CH --> C + H ), 
    (@reaction n_H * k872, CN⁻ --> CN + e ), 
    (@reaction n_H * k873, CN --> N + C ), 
    (@reaction n_H * k874, CNO --> CN + O ), 
    (@reaction n_H * k875, CO2 --> CO + O ), 
    (@reaction n_H * k876, CO --> O + C ), 
    (@reaction n_H * k877, CP --> C + P ), 
    (@reaction n_H * k878, CS --> CS⁺ + e ), 
    (@reaction n_H * k879, CS --> S + C ), 
    (@reaction n_H * k880, Cl --> Cl⁺ + e ), 
    (@reaction n_H * k881, ClO --> Cl + O ), 
    (@reaction n_H * k882, Fe --> Fe⁺ + e ), 
    (@reaction n_H * k883, H⁻ --> H + e ), 
    (@reaction n_H * k884, H2CCC --> C3H + H ), 
    (@reaction n_H * k885, H2CN --> HCN + H ), 
    (@reaction n_H * k886, H2CO --> CO + H2 ), 
    (@reaction n_H * k887, H2CS --> CS + H2 ), 
    (@reaction n_H * k888, H2O2 --> OH + OH ), 
    (@reaction n_H * k889, H2O --> OH + H ), 
    (@reaction n_H * k890, H2S2 --> HS + HS ), 
    (@reaction n_H * k891, H2S --> H2S⁺ + e ), 
    (@reaction n_H * k892, H2S --> S + H2 ), 
    (@reaction n_H * k893, H2SiO --> SiO + H2 ), 
    (@reaction n_H * k894, H --> H⁺ + e ), 
    (@reaction n_H * k895, HC2P --> CCP + H ), 
    (@reaction n_H * k896, HC3N --> CN + C2H ), 
    (@reaction n_H * k897, HC4H --> C2H + C2H ), 
    (@reaction n_H * k898, HC4H --> C4H2⁺ + e ), 
    (@reaction n_H * k899, HC4H --> C4H + H ), 
    (@reaction n_H * k900, HC5N --> C4H + CN ), 
    (@reaction n_H * k901, HC5N --> C5N + H ), 
    (@reaction n_H * k902, HC7N --> C6H + CN ), 
    (@reaction n_H * k903, HC9N --> C8H + CN ), 
    (@reaction n_H * k904, HCN --> CN + H ), 
    (@reaction n_H * k905, HCNO --> CH + NO ), 
    (@reaction n_H * k906, HCO --> CO + H ), 
    (@reaction n_H * k907, HCO --> HCO⁺ + e ), 
    (@reaction n_H * k908, HCOOCH3 --> CO2 + CH4 ), 
    (@reaction n_H * k909, HCOOCH3 --> COOCH4⁺ + e ), 
    (@reaction n_H * k910, HCOOH --> HCO + OH ), 
    (@reaction n_H * k911, HCP --> CP + H ), 
    (@reaction n_H * k912, HCS --> HCS⁺ + e ), 
    (@reaction n_H * k913, HCSi --> CH + Si ), 
    (@reaction n_H * k914, HCl --> Cl + H ), 
    (@reaction n_H * k915, HF --> H + F ), 
    (@reaction n_H * k916, HNC3 --> C2H + CN ), 
    (@reaction n_H * k917, HNC --> CN + H ), 
    (@reaction n_H * k918, HNCO --> NH + CO ), 
    (@reaction n_H * k919, HNO --> NO + H ), 
    (@reaction n_H * k920, HNSi --> SiN + H ), 
    (@reaction n_H * k921, HOCN --> OH + CN ), 
    (@reaction n_H * k922, HONC --> CN + OH ), 
    (@reaction n_H * k923, HPO --> PO + H ), 
    (@reaction n_H * k924, HS2 --> HS + S ), 
    (@reaction n_H * k925, HS --> S + H ), 
    (@reaction n_H * k926, He --> He⁺ + e ), 
    (@reaction n_H * k927, Mg --> Mg⁺ + e ), 
    (@reaction n_H * k928, N2 --> N + N ), 
    (@reaction n_H * k929, N2O --> NO + N ), 
    (@reaction n_H * k930, N --> N⁺ + e ), 
    (@reaction n_H * k931, NCCN --> CN + CN ), 
    (@reaction n_H * k932, NH2 --> NH2⁺ + e ), 
    (@reaction n_H * k933, NH2 --> NH + H ), 
    (@reaction n_H * k934, NH2CN --> NH2 + CN ), 
    (@reaction n_H * k935, NH3 --> NH2 + H ), 
    (@reaction n_H * k936, NH3 --> NH3⁺ + e ), 
    (@reaction n_H * k937, NH3 --> NH + H2 ), 
    (@reaction n_H * k938, NH --> N + H ), 
    (@reaction n_H * k939, NH --> NH⁺ + e ), 
    (@reaction n_H * k940, NO2 --> NO + O ), 
    (@reaction n_H * k941, NO --> NO⁺ + e ), 
    (@reaction n_H * k942, NO --> O + N ), 
    (@reaction n_H * k943, NS --> S + N ), 
    (@reaction n_H * k944, Na --> Na⁺ + e ), 
    (@reaction n_H * k945, O⁻ --> O + e ), 
    (@reaction n_H * k946, O2⁻ --> O2 + e ), 
    (@reaction n_H * k947, O2 --> O2⁺ + e ), 
    (@reaction n_H * k948, O2 --> O + O ), 
    (@reaction n_H * k949, O2H --> O2 + H ), 
    (@reaction n_H * k950, O --> O⁺ + e ), 
    (@reaction n_H * k951, OCN --> CN + O ), 
    (@reaction n_H * k952, OCS --> OCS⁺ + e ), 
    (@reaction n_H * k953, OCS --> S + CO ), 
    (@reaction n_H * k954, OH⁻ --> OH + e ), 
    (@reaction n_H * k955, OH --> O + H ), 
    (@reaction n_H * k956, P --> P⁺ + e ), 
    (@reaction n_H * k957, PH2 --> PH + H ), 
    (@reaction n_H * k958, PH --> P + H ), 
    (@reaction n_H * k959, PN --> P + N ), 
    (@reaction n_H * k960, PO --> P + O ), 
    (@reaction n_H * k961, S⁻ --> S + e ), 
    (@reaction n_H * k962, S2 --> S + S ), 
    (@reaction n_H * k963, S --> S⁺ + e ), 
    (@reaction n_H * k964, SO2 --> SO + O ), 
    (@reaction n_H * k965, SO --> S + O ), 
    (@reaction n_H * k966, SO --> SO⁺ + e ), 
    (@reaction n_H * k967, Si --> Si⁺ + e ), 
    (@reaction n_H * k968, SiC2 --> SiC + C ), 
    (@reaction n_H * k969, SiC2H2 --> SiC2 + H2 ), 
    (@reaction n_H * k970, SiC2H --> SiC2 + H ), 
    (@reaction n_H * k971, SiC3 --> SiC2 + C ), 
    (@reaction n_H * k972, SiC3H --> SiC3 + H ), 
    (@reaction n_H * k973, SiC4 --> SiC2 + C2 ), 
    (@reaction n_H * k974, SiC --> Si + C ), 
    (@reaction n_H * k975, SiCH2 --> SiC + H2 ), 
    (@reaction n_H * k976, SiCH3 --> SiCH2 + H ), 
    (@reaction n_H * k977, SiH2 --> SiH + H ), 
    (@reaction n_H * k978, SiH3 --> SiH2 + H ), 
    (@reaction n_H * k979, SiH4 --> SiH2 + H2 ), 
    (@reaction n_H * k980, SiH --> Si + H ), 
    (@reaction n_H * k981, SiN --> Si + N ), 
    (@reaction n_H * k982, SiNC --> Si + CN ), 
    (@reaction n_H * k983, SiO2 --> SiO + O ), 
    (@reaction n_H * k984, SiO --> Si + O ), 
    (@reaction n_H * k985, SiS --> S + Si ), 
    (@reaction n_H * k986, C10⁺ + e --> C8 + C2 ), 
    (@reaction n_H * k987, C10⁺ + e --> C9 + C ), 
    (@reaction n_H * k988, C10H⁺ + e --> C10 + H ), 
    (@reaction n_H * k989, C10H⁺ + e --> C5H + C5 ), 
    (@reaction n_H * k990, C10H⁺ + e --> C7H + C3 ), 
    (@reaction n_H * k991, C10H2⁺ + e --> C10H + H ), 
    (@reaction n_H * k992, C10H3⁺ + e --> C10H2 + H ), 
    (@reaction n_H * k993, C10H3⁺ + e --> C9H2 + C + H ), 
    (@reaction n_H * k994, C11⁺ + e --> C6 + C5 ), 
    (@reaction n_H * k995, C11⁺ + e --> C7 + C4 ), 
    (@reaction n_H * k996, C11⁺ + e --> C8 + C3 ), 
    (@reaction n_H * k997, C11⁺ + e --> C9 + C2 ), 
    (@reaction n_H * k998, C2⁺ + e --> C + C ), 
    (@reaction n_H * k999, C2H⁺ + e --> C2 + H ), 
    (@reaction n_H * k1000, C2H⁺ + e --> CH + C )] 


### Turn the Network into a system of ODEs ###
@named system = ReactionSystem(reaction_equations, t)
sys = convert(ODESystem, complete(system))
ssys = structural_simplify(sys)
prob = ODEProblem(ssys, u0, tspan, params)
sol = solve(prob, Rodas4())


### Timing ###
print("Time to convert:")
@time convert(ODESystem, complete(system))
print("Time to simplify:")
@time structural_simplify(sys)
print("Time to create the simplified problem:")
@time ODEProblem(ssys, u0, tspan, params)
print("Time to solve the simplified 1000 reaction system with Rodas4(): ")
@time solve(prob, Rodas4());

⁺
⁻


### Plotting ###
# Species number 1: C⁻
plot(sol, idxs = (0,1), lw = 3, lc = "blue", title = "Umist: Abundance of C⁻")

# Species number 2: C
plot(sol, idxs = (0,2), lw = 3, lc = "blue", title = "Umist: Abundance of C")

# Species number 3: C2
plot(sol, idxs = (0,3), lw = 3, lc = "blue", title = "Umist: Abundance of C2")

# Species number 4: e
plot(sol, idxs = (0,4), lw = 3, lc = "blue", title = "Umist: Abundance of e")

# Species number 5: CH2
plot(sol, idxs = (0,5), lw = 3, lc = "blue", title = "Umist: Abundance of CH2")

# Species number 6: C2H2
plot(sol, idxs = (0,6), lw = 3, lc = "blue", title = "Umist: Abundance of C2H2")

# Species number 7: CH
plot(sol, idxs = (0,7), lw = 3, lc = "blue", title = "Umist: Abundance of CH")

# Species number 8: C2H
plot(sol, idxs = (0,8), lw = 3, lc = "blue", title = "Umist: Abundance of C2H")

# Species number 9: CO2
plot(sol, idxs = (0,9), lw = 3, lc = "blue", title = "Umist: Abundance of CO2")

# Species number 10: CO
plot(sol, idxs = (0,10), lw = 3, lc = "blue", title = "Umist: Abundance of CO")

# Species number 11: H2O
plot(sol, idxs = (0,11), lw = 3, lc = "blue", title = "Umist: Abundance of H2O")

# Species number 12: H2CO
plot(sol, idxs = (0,12), lw = 3, lc = "blue", title = "Umist: Abundance of H2CO")

# Species number 13: N
plot(sol, idxs = (0,13), lw = 3, lc = "blue", title = "Umist: Abundance of N")

# Species number 14: CN
plot(sol, idxs = (0,14), lw = 3, lc = "blue", title = "Umist: Abundance of CN")

# Species number 15: NH
plot(sol, idxs = (0,15), lw = 3, lc = "blue", title = "Umist: Abundance of NH")

# Species number 16: HCN
plot(sol, idxs = (0,16), lw = 3, lc = "blue", title = "Umist: Abundance of HCN")

# Species number 17: O2
plot(sol, idxs = (0,17), lw = 3, lc = "blue", title = "Umist: Abundance of O2")

# Species number 18: O
plot(sol, idxs = (0,18), lw = 3, lc = "blue", title = "Umist: Abundance of O")

# Species number 19: OH
plot(sol, idxs = (0,19), lw = 3, lc = "blue", title = "Umist: Abundance of OH")

# Species number 20: HCO
plot(sol, idxs = (0,20), lw = 3, lc = "blue", title = "Umist: Abundance of HCO")

# Species number 21: C2⁻
plot(sol, idxs = (0,21), lw = 3, lc = "blue", title = "Umist: Abundance of C2⁻")

# Species number 22: C4
plot(sol, idxs = (0,22), lw = 3, lc = "blue", title = "Umist: Abundance of C4")

# Species number 23: C3
plot(sol, idxs = (0,23), lw = 3, lc = "blue", title = "Umist: Abundance of C3")

# Species number 24: C5
plot(sol, idxs = (0,24), lw = 3, lc = "blue", title = "Umist: Abundance of C5")

# Species number 25: C6
plot(sol, idxs = (0,25), lw = 3, lc = "blue", title = "Umist: Abundance of C6")

# Species number 26: C7
plot(sol, idxs = (0,26), lw = 3, lc = "blue", title = "Umist: Abundance of C7")

# Species number 27: C8
plot(sol, idxs = (0,27), lw = 3, lc = "blue", title = "Umist: Abundance of C8")

# Species number 28: C9
plot(sol, idxs = (0,28), lw = 3, lc = "blue", title = "Umist: Abundance of C9")

# Species number 29: C10
plot(sol, idxs = (0,29), lw = 3, lc = "blue", title = "Umist: Abundance of C10")

# Species number 30: C3⁻
plot(sol, idxs = (0,30), lw = 3, lc = "blue", title = "Umist: Abundance of C3⁻")

# Species number 31: C4⁻
plot(sol, idxs = (0,31), lw = 3, lc = "blue", title = "Umist: Abundance of C4⁻")

# Species number 32: C5⁻
plot(sol, idxs = (0,32), lw = 3, lc = "blue", title = "Umist: Abundance of C5⁻")

# Species number 33: C6⁻
plot(sol, idxs = (0,33), lw = 3, lc = "blue", title = "Umist: Abundance of C6⁻")

# Species number 34: C7⁻
plot(sol, idxs = (0,34), lw = 3, lc = "blue", title = "Umist: Abundance of C7⁻")

# Species number 35: C8⁻
plot(sol, idxs = (0,35), lw = 3, lc = "blue", title = "Umist: Abundance of C8⁻")

# Species number 36: C10⁻
plot(sol, idxs = (0,36), lw = 3, lc = "blue", title = "Umist: Abundance of C10⁻")

# Species number 37: C2H⁻
plot(sol, idxs = (0,37), lw = 3, lc = "blue", title = "Umist: Abundance of C2H⁻")

# Species number 38: C3H
plot(sol, idxs = (0,38), lw = 3, lc = "blue", title = "Umist: Abundance of C3H")

# Species number 39: C3H⁻
plot(sol, idxs = (0,39), lw = 3, lc = "blue", title = "Umist: Abundance of C3H⁻")

# Species number 40: C4H
plot(sol, idxs = (0,40), lw = 3, lc = "blue", title = "Umist: Abundance of C4H")

# Species number 41: C3N⁻
plot(sol, idxs = (0,41), lw = 3, lc = "blue", title = "Umist: Abundance of C3N⁻")

# Species number 42: C4N
plot(sol, idxs = (0,42), lw = 3, lc = "blue", title = "Umist: Abundance of C4N")

# Species number 43: C4H⁻
plot(sol, idxs = (0,43), lw = 3, lc = "blue", title = "Umist: Abundance of C4H⁻")

# Species number 44: C5H
plot(sol, idxs = (0,44), lw = 3, lc = "blue", title = "Umist: Abundance of C5H")

# Species number 45: C5H⁻
plot(sol, idxs = (0,45), lw = 3, lc = "blue", title = "Umist: Abundance of C5H⁻")

# Species number 46: C6H
plot(sol, idxs = (0,46), lw = 3, lc = "blue", title = "Umist: Abundance of C6H")

# Species number 47: C6H⁻
plot(sol, idxs = (0,47), lw = 3, lc = "blue", title = "Umist: Abundance of C6H⁻")

# Species number 48: C7H
plot(sol, idxs = (0,48), lw = 3, lc = "blue", title = "Umist: Abundance of C7H")

# Species number 49: C7H⁻
plot(sol, idxs = (0,49), lw = 3, lc = "blue", title = "Umist: Abundance of C7H⁻")

# Species number 50: C8H
plot(sol, idxs = (0,50), lw = 3, lc = "blue", title = "Umist: Abundance of C8H")

# Species number 51: C8H⁻
plot(sol, idxs = (0,51), lw = 3, lc = "blue", title = "Umist: Abundance of C8H⁻")

# Species number 52: C9H
plot(sol, idxs = (0,52), lw = 3, lc = "blue", title = "Umist: Abundance of C9H")

# Species number 53: C9⁻
plot(sol, idxs = (0,53), lw = 3, lc = "blue", title = "Umist: Abundance of C9⁻")

# Species number 54: C9H⁻
plot(sol, idxs = (0,54), lw = 3, lc = "blue", title = "Umist: Abundance of C9H⁻")

# Species number 55: C10H
plot(sol, idxs = (0,55), lw = 3, lc = "blue", title = "Umist: Abundance of C10H")

# Species number 56: CH⁻
plot(sol, idxs = (0,56), lw = 3, lc = "blue", title = "Umist: Abundance of CH⁻")

# Species number 57: O⁻
plot(sol, idxs = (0,57), lw = 3, lc = "blue", title = "Umist: Abundance of O⁻")

# Species number 58: OH⁻
plot(sol, idxs = (0,58), lw = 3, lc = "blue", title = "Umist: Abundance of OH⁻")

# Species number 59: S⁻
plot(sol, idxs = (0,59), lw = 3, lc = "blue", title = "Umist: Abundance of S⁻")

# Species number 60: CS
plot(sol, idxs = (0,60), lw = 3, lc = "blue", title = "Umist: Abundance of CS")

# Species number 61: CH3
plot(sol, idxs = (0,61), lw = 3, lc = "blue", title = "Umist: Abundance of CH3")

# Species number 62: CN⁻
plot(sol, idxs = (0,62), lw = 3, lc = "blue", title = "Umist: Abundance of CN⁻")

# Species number 63: CH3CN
plot(sol, idxs = (0,63), lw = 3, lc = "blue", title = "Umist: Abundance of CH3CN")

# Species number 64: CH3OH
plot(sol, idxs = (0,64), lw = 3, lc = "blue", title = "Umist: Abundance of CH3OH")

# Species number 65: HCO⁺
plot(sol, idxs = (0,65), lw = 3, lc = "blue", title = "Umist: Abundance of HCO⁺")

# Species number 66: OCS
plot(sol, idxs = (0,66), lw = 3, lc = "blue", title = "Umist: Abundance of OCS")

# Species number 67: H⁻
plot(sol, idxs = (0,67), lw = 3, lc = "blue", title = "Umist: Abundance of H⁻")

# Species number 68: CH4
plot(sol, idxs = (0,68), lw = 3, lc = "blue", title = "Umist: Abundance of CH4")

# Species number 69: H
plot(sol, idxs = (0,69), lw = 3, lc = "blue", title = "Umist: Abundance of H")

# Species number 70: H2
plot(sol, idxs = (0,70), lw = 3, lc = "blue", title = "Umist: Abundance of H2")

# Species number 71: NH2
plot(sol, idxs = (0,71), lw = 3, lc = "blue", title = "Umist: Abundance of NH2")

# Species number 72: NH3
plot(sol, idxs = (0,72), lw = 3, lc = "blue", title = "Umist: Abundance of NH3")

# Species number 73: C10H⁻
plot(sol, idxs = (0,73), lw = 3, lc = "blue", title = "Umist: Abundance of C10H⁻")

# Species number 74: C10H2
plot(sol, idxs = (0,74), lw = 3, lc = "blue", title = "Umist: Abundance of C10H2")

# Species number 75: H2CCC
plot(sol, idxs = (0,75), lw = 3, lc = "blue", title = "Umist: Abundance of H2CCC")

# Species number 76: HC3N
plot(sol, idxs = (0,76), lw = 3, lc = "blue", title = "Umist: Abundance of HC3N")

# Species number 77: HC4H
plot(sol, idxs = (0,77), lw = 3, lc = "blue", title = "Umist: Abundance of HC4H")

# Species number 78: C5H2
plot(sol, idxs = (0,78), lw = 3, lc = "blue", title = "Umist: Abundance of C5H2")

# Species number 79: C5N⁻
plot(sol, idxs = (0,79), lw = 3, lc = "blue", title = "Umist: Abundance of C5N⁻")

# Species number 80: HC5N
plot(sol, idxs = (0,80), lw = 3, lc = "blue", title = "Umist: Abundance of HC5N")

# Species number 81: C6H2
plot(sol, idxs = (0,81), lw = 3, lc = "blue", title = "Umist: Abundance of C6H2")

# Species number 82: C7H2
plot(sol, idxs = (0,82), lw = 3, lc = "blue", title = "Umist: Abundance of C7H2")

# Species number 83: C8H2
plot(sol, idxs = (0,83), lw = 3, lc = "blue", title = "Umist: Abundance of C8H2")

# Species number 84: C9H2
plot(sol, idxs = (0,84), lw = 3, lc = "blue", title = "Umist: Abundance of C9H2")

# Species number 85: HS
plot(sol, idxs = (0,85), lw = 3, lc = "blue", title = "Umist: Abundance of HS")

# Species number 86: C2N
plot(sol, idxs = (0,86), lw = 3, lc = "blue", title = "Umist: Abundance of C2N")

# Species number 87: C3N
plot(sol, idxs = (0,87), lw = 3, lc = "blue", title = "Umist: Abundance of C3N")

# Species number 88: C5N
plot(sol, idxs = (0,88), lw = 3, lc = "blue", title = "Umist: Abundance of C5N")

# Species number 89: C7N
plot(sol, idxs = (0,89), lw = 3, lc = "blue", title = "Umist: Abundance of C7N")

# Species number 90: HC7N
plot(sol, idxs = (0,90), lw = 3, lc = "blue", title = "Umist: Abundance of HC7N")

# Species number 91: C9N
plot(sol, idxs = (0,91), lw = 3, lc = "blue", title = "Umist: Abundance of C9N")

# Species number 92: HC9N
plot(sol, idxs = (0,92), lw = 3, lc = "blue", title = "Umist: Abundance of HC9N")

# Species number 93: NO
plot(sol, idxs = (0,93), lw = 3, lc = "blue", title = "Umist: Abundance of NO")

# Species number 94: NS
plot(sol, idxs = (0,94), lw = 3, lc = "blue", title = "Umist: Abundance of NS")

# Species number 95: NO2
plot(sol, idxs = (0,95), lw = 3, lc = "blue", title = "Umist: Abundance of NO2")

# Species number 96: SO2
plot(sol, idxs = (0,96), lw = 3, lc = "blue", title = "Umist: Abundance of SO2")

# Species number 97: C2O
plot(sol, idxs = (0,97), lw = 3, lc = "blue", title = "Umist: Abundance of C2O")

# Species number 98: SO
plot(sol, idxs = (0,98), lw = 3, lc = "blue", title = "Umist: Abundance of SO")

# Species number 99: H⁺
plot(sol, idxs = (0,99), lw = 3, lc = "blue", title = "Umist: Abundance of H⁺")

# Species number 100: HNC
plot(sol, idxs = (0,100), lw = 3, lc = "blue", title = "Umist: Abundance of HNC")

# Species number 101: HOC⁺
plot(sol, idxs = (0,101), lw = 3, lc = "blue", title = "Umist: Abundance of HOC⁺")

# Species number 102: O2⁻
plot(sol, idxs = (0,102), lw = 3, lc = "blue", title = "Umist: Abundance of O2⁻")

# Species number 103: C⁺
plot(sol, idxs = (0,103), lw = 3, lc = "blue", title = "Umist: Abundance of C⁺")

# Species number 104: C10H⁺
plot(sol, idxs = (0,104), lw = 3, lc = "blue", title = "Umist: Abundance of C10H⁺")

# Species number 105: C2H4
plot(sol, idxs = (0,105), lw = 3, lc = "blue", title = "Umist: Abundance of C2H4")

# Species number 106: C2H4⁺
plot(sol, idxs = (0,106), lw = 3, lc = "blue", title = "Umist: Abundance of C2H4⁺")

# Species number 107: C2H5
plot(sol, idxs = (0,107), lw = 3, lc = "blue", title = "Umist: Abundance of C2H5")

# Species number 108: C2H5⁺
plot(sol, idxs = (0,108), lw = 3, lc = "blue", title = "Umist: Abundance of C2H5⁺")

# Species number 109: C2H5OH
plot(sol, idxs = (0,109), lw = 3, lc = "blue", title = "Umist: Abundance of C2H5OH")

# Species number 110: C2H5OH⁺
plot(sol, idxs = (0,110), lw = 3, lc = "blue", title = "Umist: Abundance of C2H5OH⁺")

# Species number 111: C2O⁺
plot(sol, idxs = (0,111), lw = 3, lc = "blue", title = "Umist: Abundance of C2O⁺")

# Species number 112: C2S
plot(sol, idxs = (0,112), lw = 3, lc = "blue", title = "Umist: Abundance of C2S")

# Species number 113: C2S⁺
plot(sol, idxs = (0,113), lw = 3, lc = "blue", title = "Umist: Abundance of C2S⁺")

# Species number 114: C3O
plot(sol, idxs = (0,114), lw = 3, lc = "blue", title = "Umist: Abundance of C3O")

# Species number 115: C3O⁺
plot(sol, idxs = (0,115), lw = 3, lc = "blue", title = "Umist: Abundance of C3O⁺")

# Species number 116: C3S
plot(sol, idxs = (0,116), lw = 3, lc = "blue", title = "Umist: Abundance of C3S")

# Species number 117: C3S⁺
plot(sol, idxs = (0,117), lw = 3, lc = "blue", title = "Umist: Abundance of C3S⁺")

# Species number 118: C4H3
plot(sol, idxs = (0,118), lw = 3, lc = "blue", title = "Umist: Abundance of C4H3")

# Species number 119: C4H3⁺
plot(sol, idxs = (0,119), lw = 3, lc = "blue", title = "Umist: Abundance of C4H3⁺")

# Species number 120: C4S
plot(sol, idxs = (0,120), lw = 3, lc = "blue", title = "Umist: Abundance of C4S")

# Species number 121: C4S⁺
plot(sol, idxs = (0,121), lw = 3, lc = "blue", title = "Umist: Abundance of C4S⁺")

# Species number 122: C6H6
plot(sol, idxs = (0,122), lw = 3, lc = "blue", title = "Umist: Abundance of C6H6")

# Species number 123: C6H6⁺
plot(sol, idxs = (0,123), lw = 3, lc = "blue", title = "Umist: Abundance of C6H6⁺")

# Species number 124: CCP
plot(sol, idxs = (0,124), lw = 3, lc = "blue", title = "Umist: Abundance of CCP")

# Species number 125: CCP⁺
plot(sol, idxs = (0,125), lw = 3, lc = "blue", title = "Umist: Abundance of CCP⁺")

# Species number 126: CCl
plot(sol, idxs = (0,126), lw = 3, lc = "blue", title = "Umist: Abundance of CCl")

# Species number 127: CCl⁺
plot(sol, idxs = (0,127), lw = 3, lc = "blue", title = "Umist: Abundance of CCl⁺")

# Species number 128: CH2⁺
plot(sol, idxs = (0,128), lw = 3, lc = "blue", title = "Umist: Abundance of CH2⁺")

# Species number 129: CH2CCH2
plot(sol, idxs = (0,129), lw = 3, lc = "blue", title = "Umist: Abundance of CH2CCH2")

# Species number 130: C3H4⁺
plot(sol, idxs = (0,130), lw = 3, lc = "blue", title = "Umist: Abundance of C3H4⁺")

# Species number 131: CH2CCH
plot(sol, idxs = (0,131), lw = 3, lc = "blue", title = "Umist: Abundance of CH2CCH")

# Species number 132: CH2CCH⁺
plot(sol, idxs = (0,132), lw = 3, lc = "blue", title = "Umist: Abundance of CH2CCH⁺")

# Species number 133: CH2CN
plot(sol, idxs = (0,133), lw = 3, lc = "blue", title = "Umist: Abundance of CH2CN")

# Species number 134: CH2CN⁺
plot(sol, idxs = (0,134), lw = 3, lc = "blue", title = "Umist: Abundance of CH2CN⁺")

# Species number 135: CH2CO
plot(sol, idxs = (0,135), lw = 3, lc = "blue", title = "Umist: Abundance of CH2CO")

# Species number 136: CH2CO⁺
plot(sol, idxs = (0,136), lw = 3, lc = "blue", title = "Umist: Abundance of CH2CO⁺")

# Species number 137: CH3CCH
plot(sol, idxs = (0,137), lw = 3, lc = "blue", title = "Umist: Abundance of CH3CCH")

# Species number 138: CH3CHCH2
plot(sol, idxs = (0,138), lw = 3, lc = "blue", title = "Umist: Abundance of CH3CHCH2")

# Species number 139: C3H6⁺
plot(sol, idxs = (0,139), lw = 3, lc = "blue", title = "Umist: Abundance of C3H6⁺")

# Species number 140: CH3CHO
plot(sol, idxs = (0,140), lw = 3, lc = "blue", title = "Umist: Abundance of CH3CHO")

# Species number 141: CH3CHO⁺
plot(sol, idxs = (0,141), lw = 3, lc = "blue", title = "Umist: Abundance of CH3CHO⁺")

# Species number 142: CH3COCH3
plot(sol, idxs = (0,142), lw = 3, lc = "blue", title = "Umist: Abundance of CH3COCH3")

# Species number 143: CH3COCH3⁺
plot(sol, idxs = (0,143), lw = 3, lc = "blue", title = "Umist: Abundance of CH3COCH3⁺")

# Species number 144: CH3OCH3
plot(sol, idxs = (0,144), lw = 3, lc = "blue", title = "Umist: Abundance of CH3OCH3")

# Species number 145: CH3OCH3⁺
plot(sol, idxs = (0,145), lw = 3, lc = "blue", title = "Umist: Abundance of CH3OCH3⁺")

# Species number 146: CH⁺
plot(sol, idxs = (0,146), lw = 3, lc = "blue", title = "Umist: Abundance of CH⁺")

# Species number 147: CP
plot(sol, idxs = (0,147), lw = 3, lc = "blue", title = "Umist: Abundance of CP")

# Species number 148: CP⁺
plot(sol, idxs = (0,148), lw = 3, lc = "blue", title = "Umist: Abundance of CP⁺")

# Species number 149: ClO
plot(sol, idxs = (0,149), lw = 3, lc = "blue", title = "Umist: Abundance of ClO")

# Species number 150: ClO⁺
plot(sol, idxs = (0,150), lw = 3, lc = "blue", title = "Umist: Abundance of ClO⁺")

# Species number 151: Fe
plot(sol, idxs = (0,151), lw = 3, lc = "blue", title = "Umist: Abundance of Fe")

# Species number 152: Fe⁺
plot(sol, idxs = (0,152), lw = 3, lc = "blue", title = "Umist: Abundance of Fe⁺")

# Species number 153: H2CO⁺
plot(sol, idxs = (0,153), lw = 3, lc = "blue", title = "Umist: Abundance of H2CO⁺")

# Species number 154: H2S
plot(sol, idxs = (0,154), lw = 3, lc = "blue", title = "Umist: Abundance of H2S")

# Species number 155: H2S⁺
plot(sol, idxs = (0,155), lw = 3, lc = "blue", title = "Umist: Abundance of H2S⁺")

# Species number 156: H2SiO
plot(sol, idxs = (0,156), lw = 3, lc = "blue", title = "Umist: Abundance of H2SiO")

# Species number 157: H2SiO⁺
plot(sol, idxs = (0,157), lw = 3, lc = "blue", title = "Umist: Abundance of H2SiO⁺")

# Species number 158: C4H2⁺
plot(sol, idxs = (0,158), lw = 3, lc = "blue", title = "Umist: Abundance of C4H2⁺")

# Species number 159: HCOOCH3
plot(sol, idxs = (0,159), lw = 3, lc = "blue", title = "Umist: Abundance of HCOOCH3")

# Species number 160: COOCH4⁺
plot(sol, idxs = (0,160), lw = 3, lc = "blue", title = "Umist: Abundance of COOCH4⁺")

# Species number 161: HCP
plot(sol, idxs = (0,161), lw = 3, lc = "blue", title = "Umist: Abundance of HCP")

# Species number 162: HCP⁺
plot(sol, idxs = (0,162), lw = 3, lc = "blue", title = "Umist: Abundance of HCP⁺")

# Species number 163: HPO
plot(sol, idxs = (0,163), lw = 3, lc = "blue", title = "Umist: Abundance of HPO")

# Species number 164: HPO⁺
plot(sol, idxs = (0,164), lw = 3, lc = "blue", title = "Umist: Abundance of HPO⁺")

# Species number 165: Mg
plot(sol, idxs = (0,165), lw = 3, lc = "blue", title = "Umist: Abundance of Mg")

# Species number 166: Mg⁺
plot(sol, idxs = (0,166), lw = 3, lc = "blue", title = "Umist: Abundance of Mg⁺")

# Species number 167: NCCN
plot(sol, idxs = (0,167), lw = 3, lc = "blue", title = "Umist: Abundance of NCCN")

# Species number 168: C2N⁺
plot(sol, idxs = (0,168), lw = 3, lc = "blue", title = "Umist: Abundance of C2N⁺")

# Species number 169: CNC⁺
plot(sol, idxs = (0,169), lw = 3, lc = "blue", title = "Umist: Abundance of CNC⁺")

# Species number 170: NH3⁺
plot(sol, idxs = (0,170), lw = 3, lc = "blue", title = "Umist: Abundance of NH3⁺")

# Species number 171: NO⁺
plot(sol, idxs = (0,171), lw = 3, lc = "blue", title = "Umist: Abundance of NO⁺")

# Species number 172: NS⁺
plot(sol, idxs = (0,172), lw = 3, lc = "blue", title = "Umist: Abundance of NS⁺")

# Species number 173: Na
plot(sol, idxs = (0,173), lw = 3, lc = "blue", title = "Umist: Abundance of Na")

# Species number 174: Na⁺
plot(sol, idxs = (0,174), lw = 3, lc = "blue", title = "Umist: Abundance of Na⁺")

# Species number 175: OCS⁺
plot(sol, idxs = (0,175), lw = 3, lc = "blue", title = "Umist: Abundance of OCS⁺")

# Species number 176: P
plot(sol, idxs = (0,176), lw = 3, lc = "blue", title = "Umist: Abundance of P")

# Species number 177: P⁺
plot(sol, idxs = (0,177), lw = 3, lc = "blue", title = "Umist: Abundance of P⁺")

# Species number 178: PH
plot(sol, idxs = (0,178), lw = 3, lc = "blue", title = "Umist: Abundance of PH")

# Species number 179: PH⁺
plot(sol, idxs = (0,179), lw = 3, lc = "blue", title = "Umist: Abundance of PH⁺")

# Species number 180: PO
plot(sol, idxs = (0,180), lw = 3, lc = "blue", title = "Umist: Abundance of PO")

# Species number 181: PO⁺
plot(sol, idxs = (0,181), lw = 3, lc = "blue", title = "Umist: Abundance of PO⁺")

# Species number 182: SO⁺
plot(sol, idxs = (0,182), lw = 3, lc = "blue", title = "Umist: Abundance of SO⁺")

# Species number 183: Si
plot(sol, idxs = (0,183), lw = 3, lc = "blue", title = "Umist: Abundance of Si")

# Species number 184: Si⁺
plot(sol, idxs = (0,184), lw = 3, lc = "blue", title = "Umist: Abundance of Si⁺")

# Species number 185: SiC2
plot(sol, idxs = (0,185), lw = 3, lc = "blue", title = "Umist: Abundance of SiC2")

# Species number 186: SiC2⁺
plot(sol, idxs = (0,186), lw = 3, lc = "blue", title = "Umist: Abundance of SiC2⁺")

# Species number 187: SiC2H
plot(sol, idxs = (0,187), lw = 3, lc = "blue", title = "Umist: Abundance of SiC2H")

# Species number 188: SiC2H⁺
plot(sol, idxs = (0,188), lw = 3, lc = "blue", title = "Umist: Abundance of SiC2H⁺")

# Species number 189: SiC3
plot(sol, idxs = (0,189), lw = 3, lc = "blue", title = "Umist: Abundance of SiC3")

# Species number 190: SiC3⁺
plot(sol, idxs = (0,190), lw = 3, lc = "blue", title = "Umist: Abundance of SiC3⁺")

# Species number 191: SiC
plot(sol, idxs = (0,191), lw = 3, lc = "blue", title = "Umist: Abundance of SiC")

# Species number 192: SiC⁺
plot(sol, idxs = (0,192), lw = 3, lc = "blue", title = "Umist: Abundance of SiC⁺")

# Species number 193: SiCH2
plot(sol, idxs = (0,193), lw = 3, lc = "blue", title = "Umist: Abundance of SiCH2")

# Species number 194: SiCH2⁺
plot(sol, idxs = (0,194), lw = 3, lc = "blue", title = "Umist: Abundance of SiCH2⁺")

# Species number 195: SiCH3
plot(sol, idxs = (0,195), lw = 3, lc = "blue", title = "Umist: Abundance of SiCH3")

# Species number 196: SiCH3⁺
plot(sol, idxs = (0,196), lw = 3, lc = "blue", title = "Umist: Abundance of SiCH3⁺")

# Species number 197: SiH2
plot(sol, idxs = (0,197), lw = 3, lc = "blue", title = "Umist: Abundance of SiH2")

# Species number 198: SiH2⁺
plot(sol, idxs = (0,198), lw = 3, lc = "blue", title = "Umist: Abundance of SiH2⁺")

# Species number 199: SiH3
plot(sol, idxs = (0,199), lw = 3, lc = "blue", title = "Umist: Abundance of SiH3")

# Species number 200: SiH3⁺
plot(sol, idxs = (0,200), lw = 3, lc = "blue", title = "Umist: Abundance of SiH3⁺")

# Species number 201: SiN
plot(sol, idxs = (0,201), lw = 3, lc = "blue", title = "Umist: Abundance of SiN")

# Species number 202: SiN⁺
plot(sol, idxs = (0,202), lw = 3, lc = "blue", title = "Umist: Abundance of SiN⁺")

# Species number 203: SiS
plot(sol, idxs = (0,203), lw = 3, lc = "blue", title = "Umist: Abundance of SiS")

# Species number 204: SiS⁺
plot(sol, idxs = (0,204), lw = 3, lc = "blue", title = "Umist: Abundance of SiS⁺")

# Species number 205: C2⁺
plot(sol, idxs = (0,205), lw = 3, lc = "blue", title = "Umist: Abundance of C2⁺")

# Species number 206: S
plot(sol, idxs = (0,206), lw = 3, lc = "blue", title = "Umist: Abundance of S")

# Species number 207: S⁺
plot(sol, idxs = (0,207), lw = 3, lc = "blue", title = "Umist: Abundance of S⁺")

# Species number 208: CN⁺
plot(sol, idxs = (0,208), lw = 3, lc = "blue", title = "Umist: Abundance of CN⁺")

# Species number 209: CO⁺
plot(sol, idxs = (0,209), lw = 3, lc = "blue", title = "Umist: Abundance of CO⁺")

# Species number 210: N2⁺
plot(sol, idxs = (0,210), lw = 3, lc = "blue", title = "Umist: Abundance of N2⁺")

# Species number 211: N2
plot(sol, idxs = (0,211), lw = 3, lc = "blue", title = "Umist: Abundance of N2")

# Species number 212: O2⁺
plot(sol, idxs = (0,212), lw = 3, lc = "blue", title = "Umist: Abundance of O2⁺")

# Species number 213: C2H⁺
plot(sol, idxs = (0,213), lw = 3, lc = "blue", title = "Umist: Abundance of C2H⁺")

# Species number 214: C2H2⁺
plot(sol, idxs = (0,214), lw = 3, lc = "blue", title = "Umist: Abundance of C2H2⁺")

# Species number 215: C2H3
plot(sol, idxs = (0,215), lw = 3, lc = "blue", title = "Umist: Abundance of C2H3")

# Species number 216: C2H3⁺
plot(sol, idxs = (0,216), lw = 3, lc = "blue", title = "Umist: Abundance of C2H3⁺")

# Species number 217: C5H2⁺
plot(sol, idxs = (0,217), lw = 3, lc = "blue", title = "Umist: Abundance of C5H2⁺")

# Species number 218: C6H2⁺
plot(sol, idxs = (0,218), lw = 3, lc = "blue", title = "Umist: Abundance of C6H2⁺")

# Species number 219: C7H2⁺
plot(sol, idxs = (0,219), lw = 3, lc = "blue", title = "Umist: Abundance of C7H2⁺")

# Species number 220: C3H3⁺
plot(sol, idxs = (0,220), lw = 3, lc = "blue", title = "Umist: Abundance of C3H3⁺")

# Species number 221: CO2⁺
plot(sol, idxs = (0,221), lw = 3, lc = "blue", title = "Umist: Abundance of CO2⁺")

# Species number 222: HC3N⁺
plot(sol, idxs = (0,222), lw = 3, lc = "blue", title = "Umist: Abundance of HC3N⁺")

# Species number 223: HCN⁺
plot(sol, idxs = (0,223), lw = 3, lc = "blue", title = "Umist: Abundance of HCN⁺")

# Species number 224: C2N2⁺
plot(sol, idxs = (0,224), lw = 3, lc = "blue", title = "Umist: Abundance of C2N2⁺")

# Species number 225: C3⁺
plot(sol, idxs = (0,225), lw = 3, lc = "blue", title = "Umist: Abundance of C3⁺")

# Species number 226: C5⁺
plot(sol, idxs = (0,226), lw = 3, lc = "blue", title = "Umist: Abundance of C5⁺")

# Species number 227: CH3CH3⁺
plot(sol, idxs = (0,227), lw = 3, lc = "blue", title = "Umist: Abundance of CH3CH3⁺")

# Species number 228: CH3CH3
plot(sol, idxs = (0,228), lw = 3, lc = "blue", title = "Umist: Abundance of CH3CH3")

# Species number 229: PN⁺
plot(sol, idxs = (0,229), lw = 3, lc = "blue", title = "Umist: Abundance of PN⁺")

# Species number 230: PN
plot(sol, idxs = (0,230), lw = 3, lc = "blue", title = "Umist: Abundance of PN")

# Species number 231: H2O⁺
plot(sol, idxs = (0,231), lw = 3, lc = "blue", title = "Umist: Abundance of H2O⁺")

# Species number 232: NH2⁺
plot(sol, idxs = (0,232), lw = 3, lc = "blue", title = "Umist: Abundance of NH2⁺")

# Species number 233: O⁺
plot(sol, idxs = (0,233), lw = 3, lc = "blue", title = "Umist: Abundance of O⁺")

# Species number 234: OH⁺
plot(sol, idxs = (0,234), lw = 3, lc = "blue", title = "Umist: Abundance of OH⁺")

# Species number 235: CH3⁺
plot(sol, idxs = (0,235), lw = 3, lc = "blue", title = "Umist: Abundance of CH3⁺")

# Species number 236: CH4⁺
plot(sol, idxs = (0,236), lw = 3, lc = "blue", title = "Umist: Abundance of CH4⁺")

# Species number 237: CH3OH⁺
plot(sol, idxs = (0,237), lw = 3, lc = "blue", title = "Umist: Abundance of CH3OH⁺")

# Species number 238: N⁺
plot(sol, idxs = (0,238), lw = 3, lc = "blue", title = "Umist: Abundance of N⁺")

# Species number 239: SO2⁺
plot(sol, idxs = (0,239), lw = 3, lc = "blue", title = "Umist: Abundance of SO2⁺")

# Species number 240: CS⁺
plot(sol, idxs = (0,240), lw = 3, lc = "blue", title = "Umist: Abundance of CS⁺")

# Species number 241: Cl⁺
plot(sol, idxs = (0,241), lw = 3, lc = "blue", title = "Umist: Abundance of Cl⁺")

# Species number 242: Cl
plot(sol, idxs = (0,242), lw = 3, lc = "blue", title = "Umist: Abundance of Cl")

# Species number 243: C10⁺
plot(sol, idxs = (0,243), lw = 3, lc = "blue", title = "Umist: Abundance of C10⁺")

# Species number 244: C10H2⁺
plot(sol, idxs = (0,244), lw = 3, lc = "blue", title = "Umist: Abundance of C10H2⁺")

# Species number 245: C3H2
plot(sol, idxs = (0,245), lw = 3, lc = "blue", title = "Umist: Abundance of C3H2")

# Species number 246: C3H2⁺
plot(sol, idxs = (0,246), lw = 3, lc = "blue", title = "Umist: Abundance of C3H2⁺")

# Species number 247: C3H⁺
plot(sol, idxs = (0,247), lw = 3, lc = "blue", title = "Umist: Abundance of C3H⁺")

# Species number 248: C4⁺
plot(sol, idxs = (0,248), lw = 3, lc = "blue", title = "Umist: Abundance of C4⁺")

# Species number 249: C4H⁺
plot(sol, idxs = (0,249), lw = 3, lc = "blue", title = "Umist: Abundance of C4H⁺")

# Species number 250: C4P
plot(sol, idxs = (0,250), lw = 3, lc = "blue", title = "Umist: Abundance of C4P")

# Species number 251: C4P⁺
plot(sol, idxs = (0,251), lw = 3, lc = "blue", title = "Umist: Abundance of C4P⁺")

# Species number 252: C5H⁺
plot(sol, idxs = (0,252), lw = 3, lc = "blue", title = "Umist: Abundance of C5H⁺")

# Species number 253: C6⁺
plot(sol, idxs = (0,253), lw = 3, lc = "blue", title = "Umist: Abundance of C6⁺")

# Species number 254: C6H⁺
plot(sol, idxs = (0,254), lw = 3, lc = "blue", title = "Umist: Abundance of C6H⁺")

# Species number 255: C7⁺
plot(sol, idxs = (0,255), lw = 3, lc = "blue", title = "Umist: Abundance of C7⁺")

# Species number 256: C7H⁺
plot(sol, idxs = (0,256), lw = 3, lc = "blue", title = "Umist: Abundance of C7H⁺")

# Species number 257: C8⁺
plot(sol, idxs = (0,257), lw = 3, lc = "blue", title = "Umist: Abundance of C8⁺")

# Species number 258: C8H2⁺
plot(sol, idxs = (0,258), lw = 3, lc = "blue", title = "Umist: Abundance of C8H2⁺")

# Species number 259: C8H⁺
plot(sol, idxs = (0,259), lw = 3, lc = "blue", title = "Umist: Abundance of C8H⁺")

# Species number 260: C9⁺
plot(sol, idxs = (0,260), lw = 3, lc = "blue", title = "Umist: Abundance of C9⁺")

# Species number 261: C9H2⁺
plot(sol, idxs = (0,261), lw = 3, lc = "blue", title = "Umist: Abundance of C9H2⁺")

# Species number 262: C9H⁺
plot(sol, idxs = (0,262), lw = 3, lc = "blue", title = "Umist: Abundance of C9H⁺")

# Species number 263: CH3C4H
plot(sol, idxs = (0,263), lw = 3, lc = "blue", title = "Umist: Abundance of CH3C4H")

# Species number 264: CH3C4H⁺
plot(sol, idxs = (0,264), lw = 3, lc = "blue", title = "Umist: Abundance of CH3C4H⁺")

# Species number 265: CH3C6H
plot(sol, idxs = (0,265), lw = 3, lc = "blue", title = "Umist: Abundance of CH3C6H")

# Species number 266: C7H4⁺
plot(sol, idxs = (0,266), lw = 3, lc = "blue", title = "Umist: Abundance of C7H4⁺")

# Species number 267: CH3CN⁺
plot(sol, idxs = (0,267), lw = 3, lc = "blue", title = "Umist: Abundance of CH3CN⁺")

# Species number 268: H2CS
plot(sol, idxs = (0,268), lw = 3, lc = "blue", title = "Umist: Abundance of H2CS")

# Species number 269: H2CS⁺
plot(sol, idxs = (0,269), lw = 3, lc = "blue", title = "Umist: Abundance of H2CS⁺")

# Species number 270: H2S2
plot(sol, idxs = (0,270), lw = 3, lc = "blue", title = "Umist: Abundance of H2S2")

# Species number 271: H2S2⁺
plot(sol, idxs = (0,271), lw = 3, lc = "blue", title = "Umist: Abundance of H2S2⁺")

# Species number 272: HC2P
plot(sol, idxs = (0,272), lw = 3, lc = "blue", title = "Umist: Abundance of HC2P")

# Species number 273: HC2P⁺
plot(sol, idxs = (0,273), lw = 3, lc = "blue", title = "Umist: Abundance of HC2P⁺")

# Species number 274: HC5N⁺
plot(sol, idxs = (0,274), lw = 3, lc = "blue", title = "Umist: Abundance of HC5N⁺")

# Species number 275: HC7N⁺
plot(sol, idxs = (0,275), lw = 3, lc = "blue", title = "Umist: Abundance of HC7N⁺")

# Species number 276: HC9N⁺
plot(sol, idxs = (0,276), lw = 3, lc = "blue", title = "Umist: Abundance of HC9N⁺")

# Species number 277: HCSi
plot(sol, idxs = (0,277), lw = 3, lc = "blue", title = "Umist: Abundance of HCSi")

# Species number 278: HCSi⁺
plot(sol, idxs = (0,278), lw = 3, lc = "blue", title = "Umist: Abundance of HCSi⁺")

# Species number 279: HCl
plot(sol, idxs = (0,279), lw = 3, lc = "blue", title = "Umist: Abundance of HCl")

# Species number 280: HCl⁺
plot(sol, idxs = (0,280), lw = 3, lc = "blue", title = "Umist: Abundance of HCl⁺")

# Species number 281: HNSi
plot(sol, idxs = (0,281), lw = 3, lc = "blue", title = "Umist: Abundance of HNSi")

# Species number 282: HNSi⁺
plot(sol, idxs = (0,282), lw = 3, lc = "blue", title = "Umist: Abundance of HNSi⁺")

# Species number 283: HS2
plot(sol, idxs = (0,283), lw = 3, lc = "blue", title = "Umist: Abundance of HS2")

# Species number 284: HS2⁺
plot(sol, idxs = (0,284), lw = 3, lc = "blue", title = "Umist: Abundance of HS2⁺")

# Species number 285: HS⁺
plot(sol, idxs = (0,285), lw = 3, lc = "blue", title = "Umist: Abundance of HS⁺")

# Species number 286: N2O
plot(sol, idxs = (0,286), lw = 3, lc = "blue", title = "Umist: Abundance of N2O")

# Species number 287: N2O⁺
plot(sol, idxs = (0,287), lw = 3, lc = "blue", title = "Umist: Abundance of N2O⁺")

# Species number 288: NH⁺
plot(sol, idxs = (0,288), lw = 3, lc = "blue", title = "Umist: Abundance of NH⁺")

# Species number 289: PH2
plot(sol, idxs = (0,289), lw = 3, lc = "blue", title = "Umist: Abundance of PH2")

# Species number 290: PH2⁺
plot(sol, idxs = (0,290), lw = 3, lc = "blue", title = "Umist: Abundance of PH2⁺")

# Species number 291: S2
plot(sol, idxs = (0,291), lw = 3, lc = "blue", title = "Umist: Abundance of S2")

# Species number 292: S2⁺
plot(sol, idxs = (0,292), lw = 3, lc = "blue", title = "Umist: Abundance of S2⁺")

# Species number 293: SiC2H2
plot(sol, idxs = (0,293), lw = 3, lc = "blue", title = "Umist: Abundance of SiC2H2")

# Species number 294: SiC2H2⁺
plot(sol, idxs = (0,294), lw = 3, lc = "blue", title = "Umist: Abundance of SiC2H2⁺")

# Species number 295: SiC3H
plot(sol, idxs = (0,295), lw = 3, lc = "blue", title = "Umist: Abundance of SiC3H")

# Species number 296: SiC3H⁺
plot(sol, idxs = (0,296), lw = 3, lc = "blue", title = "Umist: Abundance of SiC3H⁺")

# Species number 297: SiC4
plot(sol, idxs = (0,297), lw = 3, lc = "blue", title = "Umist: Abundance of SiC4")

# Species number 298: SiC4⁺
plot(sol, idxs = (0,298), lw = 3, lc = "blue", title = "Umist: Abundance of SiC4⁺")

# Species number 299: SiH4
plot(sol, idxs = (0,299), lw = 3, lc = "blue", title = "Umist: Abundance of SiH4")

# Species number 300: SiH4⁺
plot(sol, idxs = (0,300), lw = 3, lc = "blue", title = "Umist: Abundance of SiH4⁺")

# Species number 301: SiH
plot(sol, idxs = (0,301), lw = 3, lc = "blue", title = "Umist: Abundance of SiH")

# Species number 302: SiH⁺
plot(sol, idxs = (0,302), lw = 3, lc = "blue", title = "Umist: Abundance of SiH⁺")

# Species number 303: SiNC
plot(sol, idxs = (0,303), lw = 3, lc = "blue", title = "Umist: Abundance of SiNC")

# Species number 304: SiNC⁺
plot(sol, idxs = (0,304), lw = 3, lc = "blue", title = "Umist: Abundance of SiNC⁺")

# Species number 305: SiO
plot(sol, idxs = (0,305), lw = 3, lc = "blue", title = "Umist: Abundance of SiO")

# Species number 306: SiO⁺
plot(sol, idxs = (0,306), lw = 3, lc = "blue", title = "Umist: Abundance of SiO⁺")

# Species number 307: H2⁺
plot(sol, idxs = (0,307), lw = 3, lc = "blue", title = "Umist: Abundance of H2⁺")

# Species number 308: F⁺
plot(sol, idxs = (0,308), lw = 3, lc = "blue", title = "Umist: Abundance of F⁺")

# Species number 309: F
plot(sol, idxs = (0,309), lw = 3, lc = "blue", title = "Umist: Abundance of F")

# Species number 310: He⁺
plot(sol, idxs = (0,310), lw = 3, lc = "blue", title = "Umist: Abundance of He⁺")

# Species number 311: He
plot(sol, idxs = (0,311), lw = 3, lc = "blue", title = "Umist: Abundance of He")

# Species number 312: C3N⁺
plot(sol, idxs = (0,312), lw = 3, lc = "blue", title = "Umist: Abundance of C3N⁺")

# Species number 313: HCNO
plot(sol, idxs = (0,313), lw = 3, lc = "blue", title = "Umist: Abundance of HCNO")

# Species number 314: HCNO⁺
plot(sol, idxs = (0,314), lw = 3, lc = "blue", title = "Umist: Abundance of HCNO⁺")

# Species number 315: HNCO
plot(sol, idxs = (0,315), lw = 3, lc = "blue", title = "Umist: Abundance of HNCO")

# Species number 316: HNCO⁺
plot(sol, idxs = (0,316), lw = 3, lc = "blue", title = "Umist: Abundance of HNCO⁺")

# Species number 317: HONC
plot(sol, idxs = (0,317), lw = 3, lc = "blue", title = "Umist: Abundance of HONC")

# Species number 318: HONC⁺
plot(sol, idxs = (0,318), lw = 3, lc = "blue", title = "Umist: Abundance of HONC⁺")

# Species number 319: HNO⁺
plot(sol, idxs = (0,319), lw = 3, lc = "blue", title = "Umist: Abundance of HNO⁺")

# Species number 320: HNO
plot(sol, idxs = (0,320), lw = 3, lc = "blue", title = "Umist: Abundance of HNO")

# Species number 321: HCOOH
plot(sol, idxs = (0,321), lw = 3, lc = "blue", title = "Umist: Abundance of HCOOH")

# Species number 322: HCOOH⁺
plot(sol, idxs = (0,322), lw = 3, lc = "blue", title = "Umist: Abundance of HCOOH⁺")

# Species number 323: NO2⁺
plot(sol, idxs = (0,323), lw = 3, lc = "blue", title = "Umist: Abundance of NO2⁺")

# Species number 324: C11
plot(sol, idxs = (0,324), lw = 3, lc = "blue", title = "Umist: Abundance of C11")

# Species number 325: C2H5CN
plot(sol, idxs = (0,325), lw = 3, lc = "blue", title = "Umist: Abundance of C2H5CN")

# Species number 326: CH2CHCNH⁺
plot(sol, idxs = (0,326), lw = 3, lc = "blue", title = "Umist: Abundance of CH2CHCNH⁺")

# Species number 327: C3P
plot(sol, idxs = (0,327), lw = 3, lc = "blue", title = "Umist: Abundance of C3P")

# Species number 328: CH2CHCCH
plot(sol, idxs = (0,328), lw = 3, lc = "blue", title = "Umist: Abundance of CH2CHCCH")

# Species number 329: CH2CHCHCH2
plot(sol, idxs = (0,329), lw = 3, lc = "blue", title = "Umist: Abundance of CH2CHCHCH2")

# Species number 330: CH2CHCN
plot(sol, idxs = (0,330), lw = 3, lc = "blue", title = "Umist: Abundance of CH2CHCN")

# Species number 331: CH2NH
plot(sol, idxs = (0,331), lw = 3, lc = "blue", title = "Umist: Abundance of CH2NH")

# Species number 332: CH2PH
plot(sol, idxs = (0,332), lw = 3, lc = "blue", title = "Umist: Abundance of CH2PH")

# Species number 333: CH3C3N
plot(sol, idxs = (0,333), lw = 3, lc = "blue", title = "Umist: Abundance of CH3C3N")

# Species number 334: CH3C5N
plot(sol, idxs = (0,334), lw = 3, lc = "blue", title = "Umist: Abundance of CH3C5N")

# Species number 335: CH3C7N
plot(sol, idxs = (0,335), lw = 3, lc = "blue", title = "Umist: Abundance of CH3C7N")

# Species number 336: CNO
plot(sol, idxs = (0,336), lw = 3, lc = "blue", title = "Umist: Abundance of CNO")

# Species number 337: H2CN
plot(sol, idxs = (0,337), lw = 3, lc = "blue", title = "Umist: Abundance of H2CN")

# Species number 338: H2O2
plot(sol, idxs = (0,338), lw = 3, lc = "blue", title = "Umist: Abundance of H2O2")

# Species number 339: HCS
plot(sol, idxs = (0,339), lw = 3, lc = "blue", title = "Umist: Abundance of HCS")

# Species number 340: HCS⁺
plot(sol, idxs = (0,340), lw = 3, lc = "blue", title = "Umist: Abundance of HCS⁺")

# Species number 341: HF
plot(sol, idxs = (0,341), lw = 3, lc = "blue", title = "Umist: Abundance of HF")

# Species number 342: HNC3
plot(sol, idxs = (0,342), lw = 3, lc = "blue", title = "Umist: Abundance of HNC3")

# Species number 343: HOCN
plot(sol, idxs = (0,343), lw = 3, lc = "blue", title = "Umist: Abundance of HOCN")

# Species number 344: NH2CN
plot(sol, idxs = (0,344), lw = 3, lc = "blue", title = "Umist: Abundance of NH2CN")

# Species number 345: O2H
plot(sol, idxs = (0,345), lw = 3, lc = "blue", title = "Umist: Abundance of O2H")

# Species number 346: OCN
plot(sol, idxs = (0,346), lw = 3, lc = "blue", title = "Umist: Abundance of OCN")

# Species number 347: SiO2
plot(sol, idxs = (0,347), lw = 3, lc = "blue", title = "Umist: Abundance of SiO2")

# Species number 348: C10H3⁺
plot(sol, idxs = (0,348), lw = 3, lc = "blue", title = "Umist: Abundance of C10H3⁺")

# Species number 349: C11⁺
plot(sol, idxs = (0,349), lw = 3, lc = "blue", title = "Umist: Abundance of C11⁺")



