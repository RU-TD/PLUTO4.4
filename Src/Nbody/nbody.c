#include "pluto.h"

double nbodySmoothingSquared(double *v, double x1, double x2, double x3)
{
    return G_SMOOTHING*G_SMOOTHING*g_isoSoundSpeed*g_isoSoundSpeed*x1*x1;
}
