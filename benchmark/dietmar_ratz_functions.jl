# Example functions from Dietmar Ratz - An Optimized Interval Slope Arithmetic and its Application

dr1(x) = (x + sin(x)) * exp(-x^2)
dr2(x) = x^4 - 10x^3 + 35x^2 - 50x + 24
dr3(x) = (log(x + 1.25) - 0.84x) ^ 2
dr4(x) = 0.02x^2 - 0.03exp(-(20(x - 0.875))^2)
dr5(x) = exp(x^2)
dr6(x) = x^4 - 12x^3 + 47x^2 - 60x - 20exp(-x)
dr7(x) = x^6 - 15x^4 + 27x^2 + 250
dr8(x) = atan(cos(tan(x)))
dr9(x) = asin(cos(acos(sin(x))))

dr_functions = [dr1, dr2, dr3, dr4, dr5, dr6, dr7, dr8, dr9]
