#include "arcball.h"
void    ArcBallSetBounds(GLfloat NewWidth, GLfloat NewHeight, ArcBall *arcball)
{
    assert((NewWidth > 1.0f) && (NewHeight > 1.0f));

    //Set adjustment factor for width/height
    arcball->AdjustWidth  = 1.0f / ((NewWidth  - 1.0f) * 0.5f);
    arcball->AdjustHeight = 1.0f / ((NewHeight - 1.0f) * 0.5f);
}


void ArcBallMapToSphere(const Point2fT* NewPt, Vector3fT* NewVec, ArcBall *arcball) 
{
    Point2fT TempPt;
    GLfloat length;

    //Copy paramter into temp point
    TempPt = *NewPt;

    //Adjust point coords and scale down to range of [-1 ... 1]
    TempPt.s.X  =        (TempPt.s.X * arcball->AdjustWidth)  - 1.0f;
    TempPt.s.Y  = 1.0f - (TempPt.s.Y * arcball->AdjustHeight);

    //Compute the square of the length of the vector to the point from the center
    length      = (TempPt.s.X * TempPt.s.X) + (TempPt.s.Y * TempPt.s.Y);

    //If the point is mapped outside of the sphere... (length > radius squared)
    if (length > 1.0f)
    {
        GLfloat norm;

        //Compute a normalizing factor (radius / sqrt(length))
        norm    = 1.0f / FuncSqrt(length);

        //Return the "normalized" vector, a point on the sphere
        NewVec->s.X = TempPt.s.X * norm;
        NewVec->s.Y = TempPt.s.Y * norm;
        NewVec->s.Z = 0.0f;
    }
    else    //Else it's on the inside
    {
        //Return a vector to a point mapped inside the sphere sqrt(radius squared - length)
        NewVec->s.X = TempPt.s.X;
        NewVec->s.Y = TempPt.s.Y;
        NewVec->s.Z = FuncSqrt(1.0f - length);
    }
}

ArcBallInit(unsigned int width,unsigned int height, ArcBall *arcball)
{
    //Clear initial values
   	GLfloat NewWidth = (GLfloat) width;
    GLfloat NewHeight = (GLfloat) height;
    arcball->StVec.s.X     =
    arcball->StVec.s.Y     = 
    arcball->StVec.s.Z     = 

    arcball->EnVec.s.X     =
    arcball->EnVec.s.Y     = 
    arcball->EnVec.s.Z     = 0.0f;

    //Set initial bounds
    ArcBallSetBounds(NewWidth, NewHeight, arcball);
}

//Mouse down
void    ArcBallClick(const Point2fT* NewPt, ArcBall *arcball)
{
    //Map the point to the sphere
    ArcBallMapToSphere(NewPt, &arcball->StVec, arcball);
}

//Mouse drag, calculate rotation
void    ArcBallDrag(const Point2fT* NewPt, Quat4fT* NewRot, ArcBall *arcball)
{
    //Map the point to the sphere
    ArcBallMapToSphere(NewPt, &arcball->EnVec, arcball);

    //Return the quaternion equivalent to the rotation
    if (NewRot)
    {
        Vector3fT  Perp;

        //Compute the vector perpendicular to the begin and end vectors
        Vector3fCross(&Perp, &arcball->StVec, &arcball->EnVec);

        //Compute the length of the perpendicular vector
        if (Vector3fLength(&Perp) > Epsilon)    //if its non-zero
        {
            //We're ok, so return the perpendicular vector as the transform after all
            NewRot->s.X = Perp.s.X;
            NewRot->s.Y = Perp.s.Y;
            NewRot->s.Z = Perp.s.Z;
            //In the quaternion values, w is cosine (theta / 2), where theta is rotation angle
            NewRot->s.W= Vector3fDot(&arcball->StVec, &arcball->EnVec);
        }
        else                                    //if its zero
        {
            //The begin and end vectors coincide, so return an identity transform
            NewRot->s.X = 
            NewRot->s.Y = 
            NewRot->s.Z = 
            NewRot->s.W = 0.0f;
        }
    }
}

