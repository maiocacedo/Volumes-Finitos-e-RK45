double x = x0;
double k1[{nfuncs}], k2[{nfuncs}], k3[{nfuncs}], k4[{nfuncs}], k5[{nfuncs}], k6[{nfuncs}];
double yn4[{nfuncs}], yn5[{nfuncs}], diff[{nfuncs}], s[{nfuncs}], yt[{nfuncs}];

for(int step = 0; step < nsteps; ++step)
{
    {
        {
            deriv_code('y', 'k1');
        }

        for (int i = 0; i < nfuncs; ++i)
            yt[i] = y[i] + 0.25 * k1[i];
        x = x0 + step * h + 0.25 * h;
        {
            deriv_code('yt', 'k2');
        }

        for (int i = 0; i < nfuncs; ++i)
            yt[i] = y[i] + (3.0 / 32.0) * k1[i] + (9.0 / 32.0) * k2[i];
        x = x0 + step * h + (3.0 / 8.0) * h;
        {
            deriv_code('yt', 'k3');
        }

        for (int i = 0; i < nfuncs; ++i)
            yt[i] = y[i] + (1932.0 / 2197.0) * k1[i] - (7200.0 / 2197.0) * k2[i] + (7296.0 / 2197.0) * k3[i];
        x = x0 + step * h + (12.0 / 13.0) * h;
        {
            deriv_code('yt', 'k4');
        }

        for (int i = 0; i < nfuncs; ++i)
            yt[i] = y[i] + (439.0 / 216.0) * k1[i] - 8.0 * k2[i] + (3680.0 / 513.0) * k3[i] - (845.0 / 4104.0) * k4[i];
        x = x0 + step * h + h;
        {
            deriv_code('yt', 'k5');
        }

        for (int i = 0; i < nfuncs; ++i)
            yt[i] = y[i] - (8.0 / 27.0) * k1[i] + 2.0 * k2[i] - (3544.0 / 2565.0) * k3[i] +
                    (1859.0 / 4104.0) * k4[i] - (11.0 / 40.0) * k5[i];
        x = x0 + step * h + 0.5 * h;
        {
            deriv_code('yt', 'k6');
        }

        for (int i = 0; i < nfuncs; ++i)
        {
            {
                yn4[i] = y[i] + (25.0 / 216.0) * k1[i] + (1408.0 / 2565.0) * k3[i] +
                         (2197.0 / 4104.0) * k4[i] - (1.0 / 5.0) * k5[i];
                yn5[i] = y[i] + (16.0 / 135.0) * k1[i] + (6656.0 / 12825.0) * k3[i] +
                         (28561.0 / 56430.0) * k4[i] - (9.0 / 50.0) * k5[i] + (2.0 / 55.0) * k6[i];
                diff[i] = fabs(yn5[i] - yn4[i]);
                s[i] = 0.84 * pow(tol * h / diff[i], 0.25);
            }
        }
        for (int i = 0; i < nfuncs; ++i)
        {
            {
                h = diff[i] > tol ? s[i] * h : h;
                y[i] = yn4[i];
            }
        }
    }
}
