// Helmholz Free Energy Equation of State
// Substance: R1234yf
// Range of validity: 220 <= T <= 410 K / p <= 30 MPA

// Reference: Markus Richter, Mark O. McLinden, and Eric W. Lemmon. Thermodynamic Properties of 2,3,3,3-Tetrafluoroprop-1-ene(R1234yf): Vapor Pressure and p-G-T Measurements and an Equation of State

function helm_r1234yf(rho, T) {
    // R1234yf constants
    var Tc = 367.85 // K
    var rhoc = 475.55 // mol.dm^-3 (475.55 kg.m^-3)
    var pc = 3382.2 // kPa

    var M = 114.04159 // g/mol
    var Rmol = 8.314472 // J/(mol K)
    var R = Rmol/(M/1000) // J/(kg K)

    var tal = Tc / T;
    var delta = rho / rhoc;

    // values for vi and ui
    var iCoef = [
        [7.549, 718.0],
        [1.537, 877.0],
        [2.030, 4465.0],
        [7.455, 1755.0]
    ];

    // coefficientes and exponents of the EoS
    // Nk | tk | dk | lk | nik | betak | gamak | epsilonk
    var kCoef = [
        [0.04592563, 1.0, 4],
        [1.546958, 0.32, 1],
        [-2.355237, 0.929, 1],
        [-0.4827835, 0.94, 2],
        [0.1758022, 0.38, 3],
        [-1.210006, 2.28, 1, 2],
        [-0.6177084, 1.76, 3, 2],
        [0.6805262, 0.97, 2, 1],
        [-0.6968555, 2.44, 2, 2],
        [-0.02695779, 1.05, 7, 1],
        [1.389966, 1.4, 1, undefined, 1.02, 1.42, 1.13, 0.712],
        [-0.4777136, 3.0, 1, undefined, 1.336, 2.31, 0.67, 0.910],
        [-0.1975184, 3.5, 3, undefined, 1.055, 0.89, 0.46, 0.677],
        [-1.147646, 1.0, 3, undefined, 5.84, 80.0, 1.28, 0.718],
        [0.0003428541, 3.5, 2, undefined, 16.2, 108.0, 1.20, 1.640]
    ];

    // Ideal gas Helmholz free energy
    // Reference state: Enthalpy = 200 kJ.kg^-1 / Entropy = 1 kJ.kg^-1.K^-1 for saturated liquid at 0 ºC

    // Ideal gas part constants:
    var a1 = -12.837928;
    var a2 = 8.042605;

    // Ideal gas part calculation
    function alpha0() {
        var alpha = 0;
        var sum = 0;
        for (var i = 0; i < 4; i++) {
            sum += iCoef[i][0] * Math.log(1 - Math.exp((-iCoef[i][1] * tal) / Tc));
        }
        alpha = a1 + a2 * tal + Math.log(delta) + 4.944 * Math.log(tal) + sum;
        return alpha;
    };

    // Ideal gas part derivatives
    function alpha0tal() {
        var alpha = 0;
        var sum = 0;

        for (var i = 0; i < 4; i++) {
            sum += iCoef[i][0] * (iCoef[i][1] / Tc) * (1 / (Math.exp(iCoef[i][1] * (tal / Tc)) - 1));
        }

        alpha = 4.944 + a2 * tal + tal * sum;
        return alpha;
    };

    function alpha0taltal() {
        var alpha = 0;
        var sum = 0;

        for (var i = 0; i < 4; i++) {
            sum += iCoef[i][0] * Math.pow((iCoef[i][1] / Tc), 2) * (Math.exp(iCoef[i][1] * (tal / Tc)) / Math.pow(Math.exp(iCoef[i][1] * (tal / Tc)), 2));
        }

        alpha = -4.944 - (tal * tal) * sum;
        return alpha;
    };

    // Residual part calculation
    function alphaR() {
        var alpha = 0;
        var sum1 = 0;
        var sum2 = 0;
        var sum3 = 0;

        for (var i = 0; i < 5; i++) {
            // CLEAR
            sum1 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]);
            // sum1 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk);
        }

        for (var i = 5; i < 10; i++) {
            // CLEAR
            sum2 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * Math.exp(-Math.pow(delta, kCoef[i][3]));
            // sum2 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * Math.exp(-Math.pow(delta, lk));
        }

        for (var i = 10; i < 15; i++) {
            // CLEAR
            sum3 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * Math.exp(-kCoef[i][4] * Math.pow((delta - kCoef[i][7]), 2) - kCoef[i][5] * Math.pow((tal - kCoef[i][6]), 2));
            // sum3 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * Math.exp(-nik * Math.pow((delta - epsilonk), 2) - betak * Math.pow((tal - gamak), 2));
        }

        alpha = sum1 + sum2 + sum3;
        return alpha;
    };

    // Residual part derivatives
    function alphaRdelta() {
        var alpha = 0;
        var sum1 = 0;
        var sum2 = 0;
        var sum3 = 0;
        
        for (var i = 0; i < 5; i++) {
            // CLEAR
            sum1 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * kCoef[i][2];
            // sum1 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * dk;
        }

        for (var i = 5; i < 10; i++) {
            // CLEAR
            sum2 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * Math.exp(-Math.pow(delta, kCoef[i][3])) * (kCoef[i][2] - kCoef[i][3] * Math.pow(delta, kCoef[i][3]));
            // sum2 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * Math.exp(-Math.pow(delta, lk)) * (dk - lk * Math.pow(delta, lk));
        }

        for (var i = 10; i < 15; i++) {
            // CLEAR
            sum3 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * Math.exp(-kCoef[i][4] * Math.pow((delta - kCoef[i][7]), 2) - kCoef[i][5] * Math.pow((tal - kCoef[i][6]), 2)) * (kCoef[i][2] - 2 * kCoef[i][4] * delta * (delta - kCoef[i][7]));
            // sum3 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * Math.exp(-nik * Math.pow((delta - epsilonk), 2) - betak * Math.pow((tal - gamak), 2)) * (dk - 2 * nik * delta * (delta - epsilonk));
        }

        alpha = sum1 + sum2 + sum3;
        return alpha;
    };

    function alphaRdeltadelta() {
        var alpha = 0;
        var sum1 = 0;
        var sum2 = 0;
        var sum3 = 0;
        
        for (var i = 0; i < 5; i++) {
            // CLEAR
            sum1 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * kCoef[i][2] * (kCoef[i][2] - 1);
            // sum1 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * dk * (dk - 1);
        }

        for (var i = 5; i < 10; i++) {
            // CLEAR
            sum2 += kCoef[i][0] * Math.pow(delta, kCoef[i][2]) * Math.pow(tal, kCoef[i][1]) * Math.exp(-Math.pow(delta, kCoef[i][3])) * ((kCoef[i][2] - kCoef[i][3] * Math.pow(delta, kCoef[i][3])) * (kCoef[i][2] - 1 - kCoef[i][3] * Math.pow(delta, kCoef[i][3])) - kCoef[i][3] * kCoef[i][3] * Math.pow(delta, kCoef[i][3]));
            // sum2 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * Math.exp(-Math.pow(delta, lk)) * ((dk - lk * Math.pow(delta, lk)) * (dk - 1 - lk * Math.pow(delta, lk)) - lk * lk * Math.pow(delta, lk));
        }

        for (var i = 10; i < 15; i++) {
            // CLEAR
            sum3 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * Math.exp(-kCoef[i][4] * Math.pow((delta - kCoef[i][7]), 2) - kCoef[i][5] * Math.pow((tal - kCoef[i][6]), 2)) * (Math.pow(kCoef[i][2] - 2 * kCoef[i][4] * delta * (delta - kCoef[i][7]), 2) - kCoef[i][2] - 2 * kCoef[i][4] * delta * delta);
        }
        alpha = sum1 + sum2 + sum3;
        return alpha;
    };

    function alphaRdeltadeltadelta() {
        var alpha = 0;
        var sum1 = 0;
        var sum2 = 0;
        var sum3 = 0;
        
        for (var i = 0; i < 5; i++) {
            // CLEAR
            sum1 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * (kCoef[i][2] * (kCoef[i][2] - 1) * (kCoef[i][2] - 2));
            // sum1 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * dk * (dk - 1) * (dk - 2);
        }

        for (var i = 5; i < 10; i++) {
            // CLEAR
            sum2 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * Math.exp(-Math.pow(delta, kCoef[i][3])) * 
            (
                kCoef[i][2] * (kCoef[i][2] - 1) * (kCoef[i][2] - 2) + kCoef[i][3] * Math.pow(delta, kCoef[i][3]) * (-2 + 6 * kCoef[i][2] - 3 * Math.pow(kCoef[i][2], 2) - 3 * kCoef[i][2] * kCoef[i][3] + 3 * kCoef[i][3] - Math.pow(kCoef[i][3], 2)) + 
                3 * Math.pow(kCoef[i][3], 2) * Math.pow(delta, 2 * kCoef[i][3]) * (kCoef[i][2] - 1 + kCoef[i][3]) - Math.pow(kCoef[i][3], 3) * Math.pow(delta, 3 * kCoef[i][3])
            );
            // sum2 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * Math.exp(-Math.pow(delta, lk)) * (dk * (dk - 1) * (dk - 2) + lk * Math.pow(delta, lk) * (-2 + 6 * dk - 3 * dk * dk - 3 * dk * lk + 3 * lk - lk * lk) + 3 * lk * lk * Math.pow(delta, 2 * lk) * (dk - 1 + lk) - Math.pow(lk, 3)*Math.pow(delta, 3*lk));
        }

        for (var i = 10; i < 15; i++) {
            // CLEAR
            sum3 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * Math.exp(-kCoef[i][4] * Math.pow((delta - kCoef[i][7]), 2) - kCoef[i][5] * Math.pow((tal - kCoef[i][6]), 2)) * 
            (
                (Math.pow(kCoef[i][2] - 2 * kCoef[i][4] * delta * (delta - kCoef[i][7]), 3)) - 3 * Math.pow(kCoef[i][2], 2) + 2 * kCoef[i][2] - 6 * kCoef[i][2] * kCoef[i][4] * delta * delta + 6 * kCoef[i][4] * delta * (delta - kCoef[i][7]) * (kCoef[i][2] + 2 * kCoef[i][4] * delta * delta)
            );
            // sum3 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * Math.exp(-nik * Math.pow((delta - epsilonk), 2)  -betak * Math.pow((tal - gamak), 2)) * (Math.pow((dk - 2 * nik * delta * (delta - epsilonk)), 3) - 3 * dk*dk + 2*dk - 6*dk*nik*delta*delta + 6*nik*delta*(delta - epsilonk)*(dk + 2*nik*delta*delta))
        }

        alpha = sum1 + sum2 + sum3;
        return alpha;
    };

    function alphaRtal() {
        var alpha = 0;
        var sum1 = 0;
        var sum2 = 0;
        var sum3 = 0;

        for (var i = 0; i < 5; i++) {
            // CLEAR
            sum1 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * kCoef[i][1];
            // sum1 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * tk;
        }

        for (var i = 5; i < 10; i++) {
            // CLEAR
            sum2 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * Math.exp(-Math.pow(delta, kCoef[i][3])) * kCoef[i][1];
            // sum2 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * Math.exp(-Math.pow(delta, lk)) * tk;
        }

        for (var i = 10; i < 15; i++) {
            // CLEAR
            sum3 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * Math.exp(-kCoef[i][4] * Math.pow((delta - kCoef[i][7]), 2) - kCoef[i][5] * Math.pow((tal - kCoef[i][6]), 2)) * (kCoef[i][1] - 2 * kCoef[i][5] * tal * (tal - kCoef[i][6]));
            // sum3 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * Math.exp(-nik * Math.pow((delta - epsilonk), 2) - betak * Math.pow((tal - gamak), 2)) * (tk - 2 * betak * tal * (tal - gamak));
        }

        alpha = sum1 + sum2 + sum3;
        return alpha;
    };

    function alphaRtaltal() {
        var alpha = 0;
        var sum1 = 0;
        var sum2 = 0;
        var sum3 = 0;

        for (var i = 0; i < 5; i++) {
            // CLEAR
            sum1 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * kCoef[i][1] * (kCoef[i][1] - 1);
            // sum1 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * tk * (tk - 1);
        }

        for (var i = 5; i < 10; i++) {
            // CLEAR
            sum2 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * Math.exp(-Math.pow(delta, kCoef[i][3])) * kCoef[i][1] * (kCoef[i][1] - 1);
            // sum2 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * Math.exp(-Math.pow(delta, lk)) * tk * (tk - 1);
        }

        for (var i = 10; i < 15; i++) {
            // CLEAR
            sum3 += kCoef[i][0] * Math.pow(delta, kCoef[i][2]) * Math.pow(tal, kCoef[i][1]) * Math.exp(-kCoef[i][4] * Math.pow((delta - kCoef[i][7]), 2) - kCoef[i][5] * Math.pow((tal - kCoef[i][6]), 2)) * (Math.pow(kCoef[i][1] - 2 * kCoef[i][5] * tal * (tal - kCoef[i][6]), 2) - kCoef[i][1] - 2 * kCoef[i][5] * tal * tal);
            // sum3 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * Math.exp(-nik * Math.pow((delta - epsilonk), 2) - betak * Math.pow((tal - gamak), 2)) * (Math.pow(tk - 2 * betak * tal * (tal - gamak), 2) - tk - 2*betak*tal*tal);
        }
        alpha = sum1 + sum2 + sum3;
        return alpha;
    };

    function alphaRtaldelta() {
        var alpha = 0;
        var sum1 = 0;
        var sum2 = 0;
        var sum3 = 0;
        
        for (var i = 0; i < 5; i++) {
            // CLEAR
            sum1 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * kCoef[i][2] * kCoef[i][1];
            // sum1 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * dk * tk;
        }

        for (var i = 5; i < 10; i++) {
            // CLEAR
            sum2 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * Math.exp(-Math.pow(delta, kCoef[i][3])) * kCoef[i][1] * (kCoef[i][2] - kCoef[i][3] * Math.pow(delta, kCoef[i][3]));
            // sum2 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * Math.exp(-Math.pow(delta, lk)) * tk * (dk - lk * Math.pow(delta, lk));
        }

        for (var i = 10; i < 15; i++) {
            // CLEAR
            sum3 += kCoef[i][0] * Math.pow(delta, kCoef[i][2]) * Math.pow(tal, kCoef[i][1]) * Math.exp(-kCoef[i][4] * Math.pow((delta - kCoef[i][7]), 2) - kCoef[i][5] * Math.pow((tal - kCoef[i][6]), 2)) * (kCoef[i][2] - 2 * kCoef[i][4] * delta * (delta - kCoef[i][7])) * (kCoef[i][1] - 2 * kCoef[i][5] * tal * (tal - kCoef[i][6]));
            // sum3 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * Math.exp(-nik * Math.pow((delta - epsilonk), 2) - betak * Math.pow((tal - gamak), 2)) * (dk - 2 * nik * delta * (delta - epsilonk)) * (tk - 2 * betak * tal * (tal - gamak));
        }

        alpha = sum1 + sum2 + sum3;
        return alpha;
    };

    function alphaRdeltataltal() {
        var alpha = 0;
        var sum1 = 0;
        var sum2 = 0;
        var sum3 = 0;

        for (var i = 0; i < 5; i++) {
            // CLEAR
            sum1 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * kCoef[i][2] * kCoef[i][1] * (kCoef[i][1] - 1);
            // sum1 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * dk * tk * (tk - 1);
        }

        for (var i = 5; i < 10; i++) {
            // CLEAR
            sum2 += kCoef[i][0] * Math.pow(tal, kCoef[i][1]) * Math.pow(delta, kCoef[i][2]) * Math.exp(-Math.pow(delta, kCoef[i][3])) * kCoef[i][1] * (kCoef[i][1] - 1) * (kCoef[i][2] - kCoef[i][3] * Math.pow(delta, kCoef[i][3]));
            // sum2 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * Math.exp(-Math.pow(delta, lk)) * (tk * (tk - 1) * (dk - lk * Math.pow(delta, lk)));
        }

        for (var i = 10; i < 15; i++) {
            // CLEAR
            sum3 += kCoef[i][0] * Math.pow(delta, kCoef[i][2]) * Math.pow(tal, kCoef[i][1]) * Math.exp(-kCoef[i][4] * Math.pow((delta - kCoef[i][7]), 2) - kCoef[i][5] * Math.pow((tal - kCoef[i][6]), 2)) * (kCoef[i][2] - 2 * kCoef[i][4] * delta * (delta - kCoef[i][7])) * (Math.pow((kCoef[i][1] - 2 * kCoef[i][5] * tal * (tal - kCoef[i][6])), 2) - kCoef[i][1] - 2 * kCoef[i][5] * tal * tal);
            // sum3 += Nk * Math.pow(delta, dk) * Math.pow(tal, tk) * Math.exp(-nik * Math.pow((delta - epsilonk), 2) - betak * Math.pow((tal - gamak), 2)) * (dk - 2 * nik * delta * (delta - epsilonk)) * (Math.pow((tk - 2 * betak * tal * (tal - gamak)), 2) - tk - 2 * betak * tal*tal);
        }

        alpha = sum1 + sum2 + sum3;
        return alpha;
    };

    return {
        cp: R * ( - (alpha0taltal() + alphaRtaltal()) + (Math.pow((1 + alphaRdelta() - alphaRtaldelta()), 2) / (1 + 2 * alphaRdelta() + alphaRdeltadelta()))),
        cv: R * (- (alpha0taltal() + alphaRtaltal())),
        h: R * T * (alpha0tal() + alphaRtal() + alphaRdelta() + 1),     // OK
        p: rho * R * T * (1 + (alphaRdelta())),                         // OK
        s: R * (alpha0tal() + alphaRtal() - alpha0() - alphaR()),       // OK
        u: R * T * (alpha0tal() + alphaRtal()),                         // OK
        v: 1/rho,
        w: (R * T / M) * Math.sqrt(1 + 2 * alphaRdelta() + alphaRdeltadelta() - (Math.pow((1 + alphaRdelta() - alphaRtaldelta()), 2) / (alpha0taltal() + alphaRtaltal())))

    };
};

// console.log(helm_r1234yf(243.231, 368.002));
console.log(helm_r1234yf(5.9, 273.15 - 30));
// helm_r1234yf(243.231, 368.002);
// helm_r1234yf(1041.014, 314.994);
// console.log(helm_r1234yf(243.231, 368.002));
// console.log(helm_r1234yf(243.231, 368.002));
// console.log(helm_r1234yf(243.231, 368.002));



