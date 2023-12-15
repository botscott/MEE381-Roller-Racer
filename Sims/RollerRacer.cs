//============================================================================
// RollerRacer.cs : Defines derived class for simulating a Roller Racer.
//       Equations of motion are derived in class notes.
//============================================================================
using Godot;
using System;

public partial class RollerRacer : Simulator
{
    // physical parameters, names are the same as that in the notes
    double m;   // mass of vehicle
    double Ig;  // moment of inertia (vertical axis) about center of mass
    double b;   // distance of com ahead of rear axle
    double c;   // distance of rear contact patch from symmetry axis
    double d;   // caster length
    double h;   // longitudinal distance between center of mass and steer axis

    double rW;  // radius of rear wheel, used for calculating rotation rate
    double rWs; // radius of steered wheel, used for calculating rotation rate

    double kPDelta;  // proportional gain for steer filter
    double kDDelta;  // derivative gain for steer filter
    double deltaDes; // desired steer angle

    double muS;      // static frict coeff, lower bound

    bool simBegun;   // indicates whether simulation has begun

    public RollerRacer() : base(11)
    {
        g = 9.81;
        muS = 0.9;
        SetInertia(25.0 /*mass*/, 0.3 /*radius of gyration*/);
        SetGeometry(1.3 /*wheel base*/, 0.6 /* cg dist from axle*/,
            0.3 /*caster dist*/, 1.0 /*wheel sep*/, 0.5 * 0.75 /*Rwheel radius*/,
            0.15 /*steered wheel radius*/);
        kPDelta = 10.0;
        kDDelta = 4.0;

        x[0] = 0.0;   // x coordinate of center of mass
        x[1] = .2;//0.0;   // xDot, time derivative of x
        x[2] = 0.0;   // z coordinate of center of mass
        x[3] = 0.0;   // zDot, time derivative of z
        x[4] = 0.0;   // psi, heading angle
        x[5] = 0.0;   // psiDot, time derivative of heading, yaw rate
        x[6] = 0.0;   // rotation angle of left rear wheel
        x[7] = 0.0;   // rotation angle of right rear wheel
        x[8] = 0.0;   // rotation angle of front steered wheel
        x[9] = 0.0;   // delta, steer angle
        x[10] = 0.0;  // deltaDot, steer rate

        SetRHSFunc(RHSFuncRRacer);
        simBegun = false;
    }

    private void RHSFuncRRacer(double[] xx, double t, double[] ff)
    {
        // give names to some state variable so code is easier to read & write
        double xDot = xx[1];
        double zDot = xx[3];
        double psi = xx[4];
        double psiDot = xx[5];
        double delta = xx[9];
        double deltaDot = xx[10];

        // calculate some trig functions here, so you only have to do it once
        double cosPsi = Math.Cos(psi);
        double sinPsi = Math.Sin(psi);
        double cosDelta = Math.Cos(delta);
        double sinDelta = Math.Sin(delta);
        double cosPsiPlusDelta = Math.Cos(psi + delta);
        double sinPsiPlusDelta = Math.Sin(psi + delta);

        // #### You will do some hefty calculations here
        // Ax-b=0
        // //x[]= {Fb        , Ff                , x..           , z..                               , psi..     , delta..   }           - B  = 0               }
        //EQ1     {sin psi   ,sin (psi + delta)  ,-m             ,0                                  ,0          ,0          }
        //EQ2     cos psi    ,cos(psi+delta)     ,0              ,-m                                 ,0          ,0          }
        //EQ3     b          ,-(h*cos(delta)-d)  ,0              ,0                                  ,-Ig        ,0          |
        //EQ7     0          ,0                  ,sin(psi)       ,cos(psi)                           ,b                                 - x. * psi. *cos(psi) - z.*psi.*sin(psi)
        //EQ10    0          ,0                  ,sin(psi+delta) ,cos(psi+delta)*cos(psi+delta)      ,0          ,0                    -z.(psi+delta)sin(psi+delta)+h*psi.*delta.*sin(delta) 
        
        double[][] A = new double[5][];

        /*
         * 
         * 0.00	0.00	-25.00	0.00	0.00	0.00	
          1.00	1.00	0.00	-25.00	0.00	0.00	
          0.60	-0.40	0.00	0.00	-2.25	0.00	
          0.00	0.00	0.00	1.00	0.60	0.00	
          0.00	0.00	0.00	1.00	0.00	-0.00
         */

        A[0] = new double[] { sinPsi, sinPsiPlusDelta, -m, 0.0, 0.0};
        A[1] = new double[] { cosPsi, cosPsiPlusDelta, 0.0, -m, 0.0 };
        A[2] = new double[] { b, -(h * cosDelta - d), 0.0, 0.0, -Ig };
        A[3] = new double[] { 0.0, 0.0, sinPsi, cosPsi, b };
        A[4] = new double[] { 0.0, 0.0, sinPsiPlusDelta, cosPsiPlusDelta * cosPsiPlusDelta, 0.0 };

        double[] b_vec = new double[] { 0.0, 0.0, 0.0, -xDot * psiDot * cosPsi + zDot * psiDot * sinPsi, -xDot*(psiDot+deltaDot)*cosPsiPlusDelta+zDot * (psi + delta)*sinPsiPlusDelta - h * psiDot * deltaDot * sinDelta };

        LinAlgEq sys2 = new LinAlgEq(5);
        sys2.A = A;
        sys2.b = b_vec;
        sys2.SolveGauss();

        // #### Right sides are zero for now. You will fix
        ff[0] = xDot; //X.
        ff[1] = sys2.sol[2]; //X..
        ff[2] = zDot;//z.
        ff[3] = sys2.sol[3];//z..
        ff[4] = psiDot;//psi.
        ff[5] = sys2.sol[4];//psi..
        ff[6] = (1.0 / rW) * (-xDot * cosPsi + zDot * sinPsi + c * psiDot); //0.0;//thetaL.
        ff[7] = (1.0 / rW) * (-xDot * cosPsi + zDot * sinPsi - c * psiDot); //0.0;//thetaR.
        ff[8] = 1.0 / rWs * (-xDot*cosPsiPlusDelta + zDot*sinPsiPlusDelta-h*psiDot*sinDelta); //0.0;//thetaF.
        ff[9] = deltaDot;
        ff[10] = -kDDelta * deltaDot - kPDelta * (delta - deltaDes);
        GD.Print(sys2.sol[2]);  //debugging
        simBegun = true;
    }

    //------------------------------------------------------------------------
    // SetInitialSpeed: Sets the initial speed of the vehicle. Must be set
    //          before simulation has begun.
    //------------------------------------------------------------------------
    public void SetInitalSpeed(double val)
    {
        if (simBegun) return;

        x[1] = val;
    }

    //------------------------------------------------------------------------
    // SetInertia: sets the two inertia properties of the vehicle. 
    //     mm: total mass in kilograms
    //     rgyr: radius of gyration in meters
    //------------------------------------------------------------------------
    public void SetInertia(double mm, double rgyr)
    {
        if (mm <= 0.1)   // check lower bound for mass
            return;     // return and not update parameters.

        if (rgyr < 0.03) // check lower bound for radius of gyration
            return;     // return and not update parameters.

        m = mm;
        Ig = m * rgyr * rgyr;
    }

    //------------------------------------------------------------------------
    // SetGeometry: Sets the geometry of the vehicle.
    //    wsb: distance between rear axle and steer axis
    //    dcg: distance from wheel axle to center of mass
    //    dcst: length of the caster
    //    wid: distance between rear wheels
    //    wRad: radius of rear wheel
    //    wRadS: radius of steered wheel
    //------------------------------------------------------------------------
    public void SetGeometry(double wsb, double dcg, double dcst, double wid,
        double wRad, double wRadS)
    {
        // check lower bounds
        if (wsb < 0.01) return;
        if (dcg <= 0.0) return;
        if (dcst < 0.0) return;
        if (wid < 0.05) return;
        if (wRad < 0.05) return;
        if (wRadS < 0.05) return;

        if (wsb - dcst < dcg) return; //cg must be btw rear axle and steer contact

        b = dcg;
        c = 0.5 * wid;
        d = dcst;
        h = wsb - dcg;

        rW = wRad;
        rWs = wRadS;
    }

    //------------------------------------------------------------------------
    // Getters/Setters
    //------------------------------------------------------------------------

    public double SteerAngleSignal
    {
        set
        {
            deltaDes = value;
        }
    }

    public double SteerAngle
    {
        get
        {
            return x[9];
        }
    }

    public double xG
    {
        get
        {
            return x[0];
        }
    }

    public double zG
    {
        get
        {
            return x[2];
        }
    }

    public double Heading
    {
        get
        {
            return x[4];
        }
    }

    public double WheelAngleL
    {
        get
        {
            return x[6];
        }
    }

    public double WheelAngleR
    {
        get
        {
            return x[7];
        }
    }

    public double WheelAngleF
    {
        get
        {
            return x[8];
        }
    }

    public double Speed
    {
        get
        {
            // ######## You have to write this part ################

            return Math.Pow(Math.Pow(x[1], 2) + Math.Pow(x[3], 2), 0.5);
        }
    }

    public double KineticEnergy
    {
        get
        {
            // ######## You have to write this part ################

            return (0.5*m*Math.Pow(Math.Pow(Math.Pow(x[1], 2) + Math.Pow(x[3], 2), 0.5),2));
        }
    }

    public double SlipRateFront
    {
        get
        {
            // ######## You have to write this part ################

            return (-1.21212121);
        }
    }

    public double SlipRateRear
    {
        get
        {
            // ######## You have to write this part ################

            return (-1.21212121);
        }
    }

    public double FontFrictionFactor
    {
        get
        {
            // ######## You have to write this part ################

            return (-1.21212121);
        }
    }
}