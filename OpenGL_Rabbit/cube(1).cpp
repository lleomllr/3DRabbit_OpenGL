/****************************************************************************************/
/*                     Projet du lapin en SYNTHESE D'IMAGE - L3 INFORMATIQUE            */
/****************************************************************************************/
/*                     Affiche a l'ecran un lapin en 3D                                */
/****************************************************************************************/

/* inclusion des fichiers d'en-tete freeglut */

#include <cstdio>
#include <cstdlib>
#include <GL/freeglut.h>
#include <jpeglib.h>
#include <iostream>
#include <cmath>

/*class Point*/
class Point{
    public :
	//coordonnées x, y et z du point
	double x;
	double y;
	double z;
};

Point pCube[8]={
    {-0.5,-0.5, 0.5},
    { 0.5, -0.5, 0.5},
    { 0.5, -0.5, -0.5},
    { -0.5, -0.5, -0.5},
    { -0.5,0.5, 0.5},
    { 0.5, 0.5, 0.5},
    { 0.5, 0.5, -0.5},
    { -0.5, 0.5, -0.5}};

//Tableau pour stocker les indices des sommets par face pour le cube
int fCube[6][4]={
  {0,3,2,1},
  {0,1,5,4},
  {1,2,6,5},
  {2,3,7,6},
  {0,4,7,3},
  {4,5,6,7}};

//initialisations d'attributs
unsigned char image[256*256*3];
const int n=6;
const double haut = 1;
char presse;
int anglex=30,angley=20,x,y,xold,yold;

float oreilAngle = 0.0f, oreilRapidite = 0.3f;
float AngleTete = 0.0, TeteRapid = 1.0;
float AngleBras = 0.0, BrasRapid = 1.0;
float AnglePieds = 0.0, PiedsRapid = 1.0;

const float PI = 3.14159265359;
Point sp[n*2];
int lum = 0;

GLuint textureID;

GLfloat RotaVertical = 0;
GLfloat RotaHorizontale = 0;
GLfloat ZoomZ = 0;
int mat = 0;

/* Prototype des fonctions */
void animationOreille();
void animationTete();
void animationBras();
void animationPieds();
void affichage();
void clavier(unsigned char touche,int x,int y);
void reshape(int x,int y);
void idle();
void mouse(int bouton,int etat,int x,int y);
void mousemotion(int x,int y);
void specialKeys(int touche, int x, int y);
void loadJpegImage(char *fichier);


//////////////////////////////////////////////////////////////////////////////////
//                                METHODES DES SOLIDES ET FONCTIONS
/////////////////////////////////////////////////////////////////////////////////

//Méthode de dessin d'une Sphère
void dessinSphere(int NM, int NP, float r){
    float coord[NP*NM][3];
    float face[(NP-1)*NM][4];

    //coordonnées des sommets :
    for(int j=0; j<NP; j++){
        for(int i=0; i<NM; i++){
            float theta = (i*2.0*PI)/(NM);
            float phi = (j*PI)/(NP-1.0) - (PI/2.0);
            coord[i+j*NM][0] = r * std::cos(theta) * std::cos(phi);
            coord[i+j*NM][1] = r* std::sin(theta)* std::cos(phi);
            coord[i+j*NM][2] = r* std::sin(phi);
        }
    }

	//indices des pts par face :
    for(int j=0; j<(NP-1); j++){
        for(int i=0; i<NM; i++){
            face[i+j*NM][0] = i + j * NM;
            face[i+j*NM][1] = ((i+1)%NM) + j* NM;
            face[i+j*NM][2] = ((i+1)%NM) + ((j+1)*NM);
                face[i+j*NM][3] = i + (j+1) * NM;
        }
    }

	//tracé de la sphère :
    for(int j=0; j<(NP-1)*(NM); j++){
        int face1 = face[j][0];
        int face2 = face[j][1];
        int face3 = face[j][2];
        int face4 = face[j][3];

        glBegin(GL_POLYGON);
        float x = j / ((NP-1)*NM*1.0);
        float y = 1;
        float z = 1;
        glVertex3f(coord[face1][0], coord[face1][1], coord[face1][2]);
        glVertex3f(coord[face2][0], coord[face2][1], coord[face2][2]);
        glVertex3f(coord[face3][0], coord[face3][1], coord[face3][2]);
        glVertex3f(coord[face4][0], coord[face4][1], coord[face4][2]);
        glEnd();
    }
}

//dessin d'une ellipse (Corps / demi ellipse pour oreilles -> nbseg/2) :
void dessinEllipse(int NM, int NP, float r){
    float coord[NP*NM][3];
    float face[(NP-1)*NM][4];

    for(int j=0; j<NP; j++){
        for(int i=0; i<NM; i++){
            float theta = (i*2.0*PI)/(NM);
            float phi = (j*PI)/(NP-1.0) - (PI/2.0);
            coord[i+j*NM][0] = r * std::cos(theta) * std::cos(phi);
            coord[i+j*NM][1] = 2*r* std::sin(theta)* std::cos(phi);
            coord[i+j*NM][2] = r* std::sin(phi);

            float s, t;
            s = (i / (float)(NM - 1));
            t = (j / (float)(NP - 1));

            glTexCoord2f(s, t);
            glVertex3f(coord[i + j * NM][0], coord[i + j * NM][1], coord[i + j * NM][2]);
        }
    }

	//indices des pts par face :
    for(int j=0; j<(NP-1); j++){
        for(int i=0; i<NM; i++){
            face[i+j*NM][0] = i + j * NM;
            face[i+j*NM][1] = ((i+1)%NM) + j* NM;
            face[i+j*NM][2] = ((i+1)%NM) + ((j+1)*NM);
                face[i+j*NM][3] = i + (j+1) * NM;
        }
    }

	//tracé de la sphère :
    for(int j=0; j<(NP-1)*(NM); j++){
        int face1 = face[j][0];
        int face2 = face[j][1];
        int face3 = face[j][2];
        int face4 = face[j][3];

        glBegin(GL_POLYGON);
        float x = j / ((NP-1)*NM*1.0);
        float y = 1;
        float z = 1;
        glTexCoord2f(coord[face1][0], coord[face1][1]);
        glVertex3f(coord[face1][0], coord[face1][1], coord[face1][2]);

        glTexCoord2f(coord[face2][0], coord[face2][1]);
        glVertex3f(coord[face2][0], coord[face2][1], coord[face2][2]);

        glTexCoord2f(coord[face3][0], coord[face3][1]);
        glVertex3f(coord[face3][0], coord[face3][1], coord[face3][2]);

        glTexCoord2f(coord[face4][0], coord[face4][1]);
        glVertex3f(coord[face4][0], coord[face4][1], coord[face4][2]);
        glEnd();
    }
}

void demieEllipse(int NM, int NP, float r){
    float coord[NP*NM][3];
    float face[(NP-1)*(NM)+NM][4];

    for(int j=0; j<NP; j++){
        for(int i=0; i<NM; i++){
            float theta = (i*2.0*PI)/(NM);
            float phi = (j*PI/2.0)/(NP-1.0) - (PI/2.0);
            coord[i+j*NM][0] = r * std::cos(theta) * std::cos(phi);
            coord[i+j*NM][1] = 2*r* std::sin(theta)* std::cos(phi);
            coord[i+j*NM][2] = r* std::sin(phi);
        }
    }

	//indices des pts par face :
    for(int j=0; j<(NP-1); j++){
        for(int i=0; i<NM; i++){
            face[i+j*NM][0] = i + j * NM;
            face[i+j*NM][1] = ((i+1)%NM) + j* NM;
            face[i+j*NM][2] = ((i+1)%NM) + ((j+1)*NM);
            face[i+j*NM][3] = i + (j+1) * NM;
        }
    }

    //faces pour le coté plat
    for (int i = 0; i < NM - 1; i++) {
        face[NM * (NP - 1) + i][0] = i;
        face[NM * (NP - 1) + i][1] = i + 1;
        face[NM * (NP - 1) + i][2] = (NM - 1) * NP + i + 1;
        face[NM * (NP - 1) + i][3] = (NM - 1) * NP + i;
    }

	//tracé de la sphère :
    for(int j=0; j<(NP-1)*(NM); j++){
        int face1 = face[j][0];
        int face2 = face[j][1];
        int face3 = face[j][2];
        int face4 = face[j][3];

        glBegin(GL_POLYGON);
        glVertex3f(coord[face1][0], coord[face1][1], coord[face1][2]);
        glVertex3f(coord[face2][0], coord[face2][1], coord[face2][2]);
        glVertex3f(coord[face3][0], coord[face3][1], coord[face3][2]);
        glVertex3f(coord[face4][0], coord[face4][1], coord[face4][2]);
        glEnd();
    }
}

//methode utile pour faire la base du nez
void demieNez(int NM, int NP, float r){
    float coord[NP*NM][3];
    float face[(NP-1)*(NM)+NM][4];

    for(int j=0; j<NP; j++){
        for(int i=0; i<NM; i++){
            float theta = (i*2.0*PI)/(NM);
            float phi = (j*PI/2.0)/(NP-1.0) - (PI/2.0);
            coord[i+j*NM][0] = r * std::cos(theta) * std::cos(phi);
            coord[i+j*NM][1] = r* std::sin(theta)* std::cos(phi);
            coord[i+j*NM][2] = 2*r* std::sin(phi);
        }
    }

	//indices des pts par face :
    for(int j=0; j<(NP-1); j++){
        for(int i=0; i<NM; i++){
            face[i+j*NM][0] = i + j * NM;
            face[i+j*NM][1] = ((i+1)%NM) + j* NM;
            face[i+j*NM][2] = ((i+1)%NM) + ((j+1)*NM);
            face[i+j*NM][3] = i + (j+1) * NM;
        }
    }

    //faces pour le coté plat
    for (int i = 0; i < NM - 1; i++) {
        face[NM * (NP - 1) + i][0] = i;
        face[NM * (NP - 1) + i][1] = i + 1;
        face[NM * (NP - 1) + i][2] = (NM - 1) * NP + i + 1;
        face[NM * (NP - 1) + i][3] = (NM - 1) * NP + i;
    }

	//tracé de la sphère :
    for(int j=0; j<(NP-1)*(NM); j++){
        int face1 = face[j][0];
        int face2 = face[j][1];
        int face3 = face[j][2];
        int face4 = face[j][3];

        glBegin(GL_POLYGON);
        glVertex3f(coord[face1][0], coord[face1][1], coord[face1][2]);
        glVertex3f(coord[face2][0], coord[face2][1], coord[face2][2]);
        glVertex3f(coord[face3][0], coord[face3][1], coord[face3][2]);
        glVertex3f(coord[face4][0], coord[face4][1], coord[face4][2]);
        glEnd();
    }
}

void dessinTete(int NM, int NP, float r){
    float coord[NP*NM][3];
    float face[(NP-1)*NM][4];

    for(int j=0; j<NP; j++){
        for(int i=0; i<NM; i++){
            float theta = (i*2.0*PI)/(NM);
            float phi = (j*PI)/(NP-1.0) - (PI/2.0);
            coord[i+j*NM][0] = 1.8*r * std::cos(theta) * std::cos(phi);
            coord[i+j*NM][1] = 2*r* std::sin(theta)* std::cos(phi);
            coord[i+j*NM][2] = 1.7*r* std::sin(phi);
        }
    }

	//indices des pts par face :
    for(int j=0; j<(NP-1); j++){
        for(int i=0; i<NM; i++){
            face[i+j*NM][0] = i + j * NM;
            face[i+j*NM][1] = ((i+1)%NM) + j* NM;
            face[i+j*NM][2] = ((i+1)%NM) + ((j+1)*NM);
                face[i+j*NM][3] = i + (j+1) * NM;
        }
    }

	//tracé de la sphère :
    for(int j=0; j<(NP-1)*(NM); j++){
        int face1 = face[j][0];
        int face2 = face[j][1];
        int face3 = face[j][2];
        int face4 = face[j][3];

        glBegin(GL_POLYGON);
        float x = j / ((NP-1)*NM*1.0);
        float y = 1;
        float z = 1;
        glVertex3f(coord[face1][0], coord[face1][1], coord[face1][2]);
        glVertex3f(coord[face2][0], coord[face2][1], coord[face2][2]);
        glVertex3f(coord[face3][0], coord[face3][1], coord[face3][2]);
        glVertex3f(coord[face4][0], coord[face4][1], coord[face4][2]);
        glEnd();
    }
}

void dessinOreille(int NM, int NP, float r){
    float coord[NP*NM][3];
    float face[(NP-1)*NM][4];

    for(int j=0; j<NP; j++){
        for(int i=0; i<NM; i++){
            float theta = (i*2.0*PI)/(NM);
            float phi = (j*PI)/(NP-1.0) - (PI/2.0);
            coord[i+j*NM][0] = 0.9*r * std::cos(theta) * std::cos(phi);
            coord[i+j*NM][1] = 5*r* std::sin(theta)* std::cos(phi);
            coord[i+j*NM][2] = 0.3*r* std::sin(phi);
        }
    }

	//indices des pts par face :
    for(int j=0; j<(NP-1); j++){
        for(int i=0; i<NM; i++){
            face[i+j*NM][0] = i + j * NM;
            face[i+j*NM][1] = ((i+1)%NM) + j* NM;
            face[i+j*NM][2] = ((i+1)%NM) + ((j+1)*NM);
                face[i+j*NM][3] = i + (j+1) * NM;
        }
    }

	//tracé de la sphère :
    for(int j=0; j<(NP-1)*(NM); j++){
        int face1 = face[j][0];
        int face2 = face[j][1];
        int face3 = face[j][2];
        int face4 = face[j][3];

        glBegin(GL_POLYGON);
        float x = j / ((NP-1)*NM*1.0);
        float y = 1;
        float z = 1;
        glVertex3f(coord[face1][0], coord[face1][1], coord[face1][2]);
        glVertex3f(coord[face2][0], coord[face2][1], coord[face2][2]);
        glVertex3f(coord[face3][0], coord[face3][1], coord[face3][2]);
        glVertex3f(coord[face4][0], coord[face4][1], coord[face4][2]);
        glEnd();
    }
}

//methode utile pour modéliser le museau du lapin
void dessinNez(int NM, int NP, float r){
    float coord[NP*NM][3];
    float face[(NP-1)*NM][4];

    for(int j=0; j<NP; j++){
        for(int i=0; i<NM; i++){
            float theta = (i*2.0*PI)/(NM);
            float phi = (j*PI)/(NP-1.0) - (PI/2.0);
            coord[i+j*NM][0] = 2*r * std::cos(theta) * std::cos(phi);
            coord[i+j*NM][1] = 1.5*r* std::sin(theta)* std::cos(phi);
            coord[i+j*NM][2] = 1.7*r* std::sin(phi);
        }
    }

	//indices des pts par face :
    for(int j=0; j<(NP-1); j++){
        for(int i=0; i<NM; i++){
            face[i+j*NM][0] = i + j * NM;
            face[i+j*NM][1] = ((i+1)%NM) + j* NM;
            face[i+j*NM][2] = ((i+1)%NM) + ((j+1)*NM);
                face[i+j*NM][3] = i + (j+1) * NM;
        }
    }

	//tracé de la sphère :
    for(int j=0; j<(NP-1)*(NM); j++){
        int face1 = face[j][0];
        int face2 = face[j][1];
        int face3 = face[j][2];
        int face4 = face[j][3];

        glBegin(GL_POLYGON);
        float x = j / ((NP-1)*NM*1.0);
        float y = 1;
        float z = 1;
        glVertex3f(coord[face1][0], coord[face1][1], coord[face1][2]);
        glVertex3f(coord[face2][0], coord[face2][1], coord[face2][2]);
        glVertex3f(coord[face3][0], coord[face3][1], coord[face3][2]);
        glVertex3f(coord[face4][0], coord[face4][1], coord[face4][2]);
        glEnd();
    }
}

void dessinOeil(int NM, int NP, float r){
    float coord[NP*NM][3];
    float face[(NP-1)*NM][4];

    for(int j=0; j<NP; j++){
        for(int i=0; i<NM; i++){
            float theta = (i*2.0*PI)/(NM);
            float phi = (j*PI)/(NP-1.0) - (PI/2.0);
            coord[i+j*NM][0] = 1.5*r * std::cos(theta) * std::cos(phi);
            coord[i+j*NM][1] = r* std::sin(theta)* std::cos(phi);
            coord[i+j*NM][2] =r* std::sin(phi);
        }
    }

	//indices des pts par face :
    for(int j=0; j<(NP-1); j++){
        for(int i=0; i<NM; i++){
            face[i+j*NM][0] = i + j * NM;
            face[i+j*NM][1] = ((i+1)%NM) + j* NM;
            face[i+j*NM][2] = ((i+1)%NM) + ((j+1)*NM);
                face[i+j*NM][3] = i + (j+1) * NM;
        }
    }

	//tracé de la sphère :
    for(int j=0; j<(NP-1)*(NM); j++){
        int face1 = face[j][0];
        int face2 = face[j][1];
        int face3 = face[j][2];
        int face4 = face[j][3];

        glBegin(GL_POLYGON);
        float x = j / ((NP-1)*NM*1.0);
        float y = 1;
        float z = 1;
        glVertex3f(coord[face1][0], coord[face1][1], coord[face1][2]);
        glVertex3f(coord[face2][0], coord[face2][1], coord[face2][2]);
        glVertex3f(coord[face3][0], coord[face3][1], coord[face3][2]);
        glVertex3f(coord[face4][0], coord[face4][1], coord[face4][2]);
        glEnd();
    }
}


//dent (petit cylindre ~détail) :
void Cylindre(float x, float y, float z)
{
	for(int i=0; i<n; i++){

            float px =x*cos((2*i*M_PI)/n);
            float pz =-z*sin((2*i*M_PI)/n);
    		Point p1;
    		p1.x = px;
    		p1.y = y/2;
    		p1.z = pz;
    		Point p2;
    		p2.x = px;
    		p2.y = -y/2;
    		p2.z=pz;
    		sp[i]=p1;
    		sp[i+n]=p2;
	}
	for(int i=0; i<n;i++){
    		glBegin(GL_POLYGON);
                glTexCoord2f(0,0);
                glVertex3f(sp[i].x,sp[i].y,sp[i].z);
                glTexCoord2f(0.5,0);

                glVertex3f(sp[(i+1)%n].x,sp[(i+1)%n].y,sp[(i+1)%n].z);
                glTexCoord2f(0.5,0.5);
                glVertex3f(sp[n+(i+1)%n].x,sp[n+(i+1)%n].y,sp[n+(i+1)%n].z);
                glTexCoord2f(0,0.5);

                glVertex3f(sp[i+n].x,sp[i+n].y,sp[i+n].z);
    		glEnd();
	}

	glBegin(GL_POLYGON);
	for(int i=0; i<n;i++){
    		glVertex3f(sp[i].x,sp[i].y,sp[i].z);
	}
	glEnd();
	glBegin(GL_POLYGON);
	for(int i=0; i<n;i++){
    		glVertex3f(sp[i+n].x,sp[i+n].y,sp[i+n].z);
	}
	glEnd();

}

//Cou du lapin sur lequel on applique une texture à partir des équations paramétriques du cylindre
void Cou()
{
    glBindTexture(GL_TEXTURE_2D, textureID);
	for(int i=0; i<n; i++){

    		double x =0.15*cos((2*i*M_PI)/n);
    		double z =-0.1*sin((2*i*M_PI)/n);
    		Point p1;
    		p1.x = x;
    		p1.y = haut/2;
    		p1.z = z;
    		Point p2;
    		p2.x = x;
    		p2.y = -haut/2;
    		p2.z=z;
    		sp[i]=p1;
    		sp[i+n]=p2;
	}
	for(int i=0; i<n;i++){
        glBegin(GL_POLYGON);

        glTexCoord2f(0, 0);
        glVertex3f(sp[i].x, sp[i].y, sp[i].z);

        glTexCoord2f(1, 0);
        glVertex3f(sp[(i + 1) % n].x, sp[(i + 1) % n].y, sp[(i + 1) % n].z);

        glTexCoord2f(1, 1);
        glVertex3f(sp[n + (i + 1) % n].x, sp[n + (i + 1) % n].y, sp[n + (i + 1) % n].z);

        glTexCoord2f(0,1);
        glVertex3f(sp[i + n].x, sp[i + n].y, sp[i + n].z);

        glEnd();
	}

	glBegin(GL_POLYGON);
	for(int i=0; i<n;i++){

    		glVertex3f(sp[i].x,sp[i].y,sp[i].z);
	}
	glEnd();
	glBegin(GL_POLYGON);
	for(int i=0; i<n;i++){

    		glVertex3f(sp[i+n].x,sp[i+n].y,sp[i+n].z);
	}
	glEnd();

	glColor3f(0,0,0);


}

//Bras (rotule dans Startbras + rotule dans finBras + Main + Doigts):
void dessinBras(){
    glPushMatrix();
    glTranslatef(0.0f,1.3f,0.0f);
    glColor3f(0.9f, 0.9f, 0.9f);
    glRotatef(90,-1,0,0);
    glutSolidCylinder(0.12,0.8,20,20);
    glPopMatrix();
}

void dessinRotule(){
    glPushMatrix();
    glTranslatef(0.0f,0.9f,0.0f);
    glutSolidSphere(0.12,20,20);
    glPopMatrix();
}

//methode qui dessine le bras gauche
void brasG()
{
    glPushMatrix();
        glTranslatef(-0.3, -0.35, 0);
       glColor3f(0.9f, 0.9f, 0.9f);
        dessinRotule();
    glPopMatrix();


    glPushMatrix();
        glTranslatef(1, 1, 0);
        glRotatef(109, 0, 0, 0.05);
        glColor3f(0.9f, 0.9f, 0.9f);
        dessinBras();
    glPopMatrix();

    glPushMatrix();
        glTranslatef(-1, -0.6, 0);
        glColor3f(0.95f, 0.95f, 0.95f);
        dessinRotule();
    glPopMatrix();

    //main gauche
    glPushMatrix();
        glTranslatef(-1.1, 0.4, 0);
        glRotatef(90, 0, 0, 1);
        dessinEllipse(20,20,0.05);
    glPopMatrix();
    glPushMatrix();
        glTranslatef(-1, 0.18, 0);
        glRotatef(140, 0, 0, 1);
        dessinEllipse(20,20,0.05);
    glPopMatrix();
    glPushMatrix();
        glTranslatef(-1.1, 0.265, 0);
        glRotatef(110, 0, 0, 1);
        dessinEllipse(20,20,0.05);
    glPopMatrix();

}

//methode qui dessine le bras droit
void brasD()
{
    glPushMatrix();
        glTranslatef(0.3, -0.35, 0);
        glColor3f(0.9f, 0.9f, 0.9f);
        dessinRotule();
    glPopMatrix();

    glPushMatrix();
        glTranslatef(-1, 1, 0);
        glRotatef(-109, 0, 0, 0.05);
        glColor3f(0.9f, 0.9f, 0.9f);
        dessinBras();
    glPopMatrix();
    glPushMatrix();
        glTranslatef(1, -0.6, 0);
        glColor3f(0.95f, 0.95f, 0.95f);
        dessinRotule();
    glPopMatrix();

    //main droite
    glPushMatrix();
        glTranslatef(1.1, 0.4, 0);
        glRotatef(-90, 0, 0, 1);
        dessinEllipse(20,20,0.05);
    glPopMatrix();
    glPushMatrix();
        glTranslatef(1, 0.18, 0);
        glRotatef(-140, 0, 0, 1);
        dessinEllipse(20,20,0.05);
    glPopMatrix();
    glPushMatrix();
        glTranslatef(1.1, 0.265, 0);
        glRotatef(-110, 0, 0, 1);
        dessinEllipse(20,20,0.05);
    glPopMatrix();
}

//base du lapin
void base()
{
     // Dessin du cube
  for (int i=0;i<6;i++)
    {
      glBegin(GL_POLYGON);
      for (int j=0;j<4;j++){
          glVertex3f(pCube[fCube[i][j]].x,pCube[fCube[i][j]].y,pCube[fCube[i][j]].z);
      }
      glEnd();
    }
}

//méthode qui modélise l'entiereté de la tête (corps de la tête, museau, nez, yeux, dents et oreilles)
void tete()
{
    //tete
    glPushMatrix();
        dessinTete(50,50,0.5);
    glPopMatrix();

    //joues
    glPushMatrix();
        glTranslatef(-0.3,-0.25,0.35);
        glColor3f(0.95f, 0.95f, 0.95f);
        dessinSphere(12,12,0.55);
    glPopMatrix();
    glPushMatrix();
        glTranslatef(0.3,-0.25,0.35);
        glColor3f(0.95f, 0.95f, 0.95f);
        dessinSphere(12,12,0.55);
    glPopMatrix();

    //museau
    glPushMatrix();
        glTranslatef(-0.17,-0.15,0.75);
        glRotatef(35, 0, 0, 0.5);
        glColor3f(0.9f, 0.9f, 0.9f);

        dessinNez(20,20,0.15);
    glPopMatrix();
    glPushMatrix();
        glTranslatef(0.17,-0.15,0.75);
        glRotatef(-35, 0, 0, 0.5);
        glColor3f(0.9f, 0.9f, 0.9f);

        dessinNez(20,20,0.15);
    glPopMatrix();

    //oreilles
    glPushMatrix();
        glTranslatef(0.5,0.75,0);
        glRotatef(-oreilAngle,0,0,1);
        glTranslatef(0,1.25,0);
        dessinOreille(20,20,0.3);

    glPopMatrix();
    glPushMatrix();
        glTranslatef(-0.5,0.75,0);
        glRotatef(oreilAngle,0,0,1);
        glTranslatef(0,1.25,0);
        dessinOreille(20,20,0.3);
    glPopMatrix();

    //yeux
    glPushMatrix();
        glTranslatef(-0.3,0.45,0.45);
        glRotatef(-90,0,0,1);
        glRotatef(-30,1,0,1);
        glRotatef(-40, 0, 1, 1);
        glScalef(1.08,1.08,1.08);
        glColor3f(0, 0, 0);
        dessinOeil(20,20,0.25);
    glPopMatrix();
    glPushMatrix();
        glTranslatef(0.3,0.45,0.45);
        glRotatef(-60,0,0,1);
        glRotatef(-50,1,0,0);
        glRotatef(-40, 0, 1, 1);
        glScalef(1.08,1.08,1.08);
        glColor3f(0, 0, 0);
        dessinOeil(20,20,0.25);
    glPopMatrix();

    //nez
    glPushMatrix();
        glTranslatef(0,0.09,0.95);
        glRotatef(270, 1, 0, 0);
        glColor3f(0,0,0);
        demieNez(20,20,0.1);
    glPopMatrix();
    glPushMatrix();
        glTranslatef(0,0.08,0.95);
        glColor3f(0,0,0);
        dessinSphere(20,20,0.1);
    glPopMatrix();

    //dent
    glPushMatrix();
        glTranslatef(-0.07,-0.3,0.9);
        glColor3f(1,1,1);
        glRotatef(90, 1, 0, 0);
        Cylindre(0.1,0.08,0.3);
    glPopMatrix();
    glPushMatrix();
        glTranslatef(0.07,-0.3,0.9);
        glColor3f(1,1,1);
        glRotatef(90, 1, 0, 0);
        Cylindre(0.1,0.08,0.3);
    glPopMatrix();

}

//methode qui gère deux types de lumières différentes
void lumiere(){
    glEnable(GL_LIGHTING);
	if(lum == 1){
		//coord homogènes : position
    		GLfloat position_source0[] = {5.0, 5.0, 5.0, 1.0};
    		//direction source à distance infinie
    		GLfloat direction_source0[] = {1.0, 2.0, 3.0, 0.0};
    		GLfloat dif_0[] = {1.0, 0.0, 0.0, 1.0};//composante diffuse rouge
    		GLfloat amb_0[] = {1.0, 0.0, 0.0, 1.0}; //composante ambiante rouge
    		GLfloat spec_0[] = {1.0, 1.0, 1.0, 1.0}; //composante spéculaire blanche
    		//spécification des propriétés
    		glLightfv(GL_LIGHT0, GL_POSITION,position_source0);
    		glLightfv(GL_LIGHT0, GL_AMBIENT, amb_0);
    		glLightfv(GL_LIGHT0, GL_DIFFUSE, dif_0);
    		glLightfv(GL_LIGHT0, GL_SPECULAR, spec_0);
    		//activation de la source GL_LIGHT0
    		glEnable(GL_LIGHT0);
	}
	if(lum == 2){
		//coord homogènes : position
    		GLfloat position_source0[] = {5.0, 5.0, 5.0, 1.0};
    		//direction source à distance infinie
    		GLfloat direction_source0[] = {0.0, 1.0, 0.0, 0.0};
    		GLfloat dif_0[] = {1.0, 1.0, 1.0, 1.0};//composante diffuse blanche
    		GLfloat amb_0[] = {0.2, 0.2, 0.2, 1.0}; //composante ambiante faible
    		GLfloat spec_0[] = {1.0, 1.0, 1.0, 1.0}; //composante spéculaire blanche
    		//spécification des propriétés
    		glLightfv(GL_LIGHT0, GL_POSITION,position_source0);
    		glLightfv(GL_LIGHT0, GL_AMBIENT, amb_0);
    		glLightfv(GL_LIGHT0, GL_DIFFUSE, dif_0);
    		glLightfv(GL_LIGHT0, GL_SPECULAR, spec_0);
    		//activation de la source GL_LIGHT0
    		glEnable(GL_LIGHT0);
	}
	if(lum == 3){
    		GLfloat position_source0[] = {5.0, 5.0, 5.0, 1.0};
    		GLfloat direction_source0[] = {0.0, 1.0, 0.0, 0.0};
    		GLfloat dif_0[] = {1.0, 1.0, 1.0, 1.0};
    		GLfloat amb_0[] = {1.0, 1.0, 1.0, 1.0};//composante ambiante élevé
    		GLfloat spec_0[] = {1.0, 1.0, 1.0, 1.0};
    		glLightfv(GL_LIGHT0, GL_POSITION,position_source0);
    		glLightfv(GL_LIGHT0, GL_AMBIENT, amb_0);
    		glLightfv(GL_LIGHT0, GL_DIFFUSE, dif_0);
    		glLightfv(GL_LIGHT0, GL_SPECULAR, spec_0);
    		glEnable(GL_LIGHT0);
	}
	if(lum == 0) glDisable(GL_LIGHTING);
}

//méthode pour appliquer des matériaux à notre lapin
void materiaux()
{
    if(mat == 1){
    GLfloat mat_ambient[] = {0.2, 0.2, 0.2, 1.0};
    GLfloat mat_diffuse[] = {0.8, 0.8, 0.8, 1.0};
    GLfloat mat_specular[] = {0, 0, 0, 1.0};
    GLfloat mat_shininess[] = {0.0};

    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    }
    if(mat == 2){
    GLfloat mat_ambient[] = {0.0, 0.0, 0.0, 1.0};
    GLfloat mat_diffuse[] = {0.5, 0.5, 0.5, 1.0};
    GLfloat mat_specular[] = {0.8, 0.8, 0.8, 1.0};
    GLfloat mat_shininess[] = {100.0};

    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    }
    if(mat == 0){
    GLfloat mat_ambient[] = {0.4, 0.2, 0.1, 1.0};
    GLfloat mat_diffuse[] = {0.4, 0.2, 0.1, 1.0};
    GLfloat mat_specular[] = {0.0, 0.0, 0.0, 1.0};
    GLfloat mat_shininess[] = {10.0};

    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    }
}


//////////////////////////////////////////////////////////////////////////////////
//                                 AFFICHAGE
//////////////////////////////////////////////////////////////////////////////////

int main(int argc,char **argv)
{
    /* Chargement de la texture */
    loadJpegImage("./texture_peau.jpg");
    glGenTextures(1, &textureID);

  /* initialisation de glut et creation
     de la fenetre */
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowPosition(200,200);
  glutInitWindowSize(1000,750);
  glutCreateWindow("Lapin");

  /* Initialisation d'OpenGL */
  glClearColor(0.0,0.0,0.0,0.0);
  glColor3f(1.0,1.0,1.0);
  glPointSize(2.0);
  glEnable(GL_DEPTH_TEST);

  /* Parametrage du placage de textures */
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
  glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,256,256,0,
	       GL_RGB,GL_UNSIGNED_BYTE,image);
  glEnable(GL_TEXTURE_2D);

  /* enregistrement des fonctions de rappel */
  glutIdleFunc(animationOreille);
  /*glutIdleFunc(animationBras);*/
  glutDisplayFunc(affichage);
  glutKeyboardFunc(clavier);
  glutReshapeFunc(reshape);
  glutMouseFunc(mouse);
  glutMotionFunc(mousemotion);
  glutSpecialFunc(specialKeys);

  /* Entree dans la boucle principale glut */
  glutMainLoop();
  return 0;
}

void affichage()
{
    lumiere();
    materiaux();
    int val = 3;

  /* effacement de l'image avec la couleur de fond */
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glShadeModel(GL_SMOOTH);

  /*glMatrixMode(GL_PROJECTION);*/
  glLoadIdentity();
  glOrtho(-val-ZoomZ, val+ZoomZ, -val-ZoomZ, val+ZoomZ, -val-ZoomZ, val+ZoomZ);
  /*glRotatef(angley,1.0,0.0,0.0);
  glRotatef(anglex,0.0,1.0,0.0);*/

  glRotatef(RotaVertical,1,0,0);
  glRotatef(RotaHorizontale,0,1,0);

    //Repère
    //axe x en rouge
    /*glBegin(GL_LINES);
        glColor3f(1.0,0.0,0.0);
    	glVertex3f(0, 0,0.0);
    	glVertex3f(1, 0,0.0);
    glEnd();
    //axe des y en vert
    glBegin(GL_LINES);
    	glColor3f(0.0,1.0,0.0);
    	glVertex3f(0, 0,0.0);
    	glVertex3f(0, 1,0.0);
    glEnd();
    //axe des z en bleu
    glBegin(GL_LINES);
    	glColor3f(0.0,0.0,1.0);
    	glVertex3f(0, 0,0.0);
    	glVertex3f(0, 0,1.0);
    glEnd();*/

  //corps
  glPushMatrix();
        glColor3f(0.9f, 0.9f, 0.85f);
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, textureID);
        dessinEllipse(20, 20, 0.5);
  glPopMatrix();
  glDisable(GL_TEXTURE_2D);

  //membre du corps
    glPushMatrix();
        brasG();
    glPopMatrix();
    glPushMatrix();
        glRotatef(AngleBras, 0, 0, 1);
        brasD();
    glPopMatrix();

    //jambes
    glPushMatrix();
        glTranslatef(0.2, -0.7, 0);
        glRotatef(-90, 0, 1, 0);
        glRotatef(-30, 0, 0, 0.5);
        glColor3f(1, 1, 1);
        dessinEllipse(20,20,0.3);
    glPopMatrix();
    glPushMatrix();
        glTranslatef(-0.2, -0.7, 0);
        glRotatef(-90, 0, 1, 0);
        glRotatef(-30, 0, 0, 0.5);
        glColor3f(1, 1, 1);
        dessinEllipse(20,20,0.3);
    glPopMatrix();

    //pieds
    glPushMatrix();
        glTranslatef(0.2,-1.3,0);
        glRotatef(-90,1, 0, 0);
        glRotatef(-180,0, 1, 0);
        glRotatef(-30,0, 0, 0.2);
        glRotatef(AnglePieds, 0, 0, 1);
        glColor3f(0.9f, 0.9f, 0.9f);
        demieEllipse(20,20,0.21);
    glPopMatrix();
    glPushMatrix();
        glTranslatef(-0.2,-1.3,0);
        glRotatef(-90,1, 0, 0);
        glRotatef(-180,0, 1, 0);
        glRotatef(30,0, 0, 0.2);
        glRotatef(-AnglePieds, 0, 0, 1);
        glColor3f(0.9f, 0.9f, 0.9f);

        demieEllipse(20,20,0.21);
    glPopMatrix();

    //queue
    glPushMatrix();
        glTranslatef(0, -0.5, -0.4);
        dessinSphere(20,20,0.2);
    glPopMatrix();

    //torse
    glPushMatrix();
        glTranslatef(0, 0.13, 0.17);
        glRotatef(-180, 0.5, 0, 0);
        glColor3f(1.0f, 0.8f, 0.8f);

        demieEllipse(30,30,0.37);
    glPopMatrix();

    //tete et mise à l'echelle
    glPushMatrix();
        glTranslatef(0,1.67,0.27);
        glScalef(0.65,0.65,0.65);
        glRotatef(AngleTete, 0, 1, 0);
        glColor3f(0.9f, 0.9f, 0.85f);
        tete();
    glPopMatrix();

    //cou
    glPushMatrix();
        glTranslatef(0,1,0);
        glEnable(GL_TEXTURE_2D);
        Cou();
        glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    //base
    glPushMatrix();
        glColor3f(0,1,0);
        glTranslatef(0,-1.34,0);
        glScalef(4,0.1,3);
        base();
    glPopMatrix();

    glFlush();
  //On echange les buffers
  glutSwapBuffers();
}


//////////////////////////////////////////////////////////////////////////////////
//                                   FONCTIONNALITES DU CLAVIER ET TEXTURES
/////////////////////////////////////////////////////////////////////////////////

void clavier(unsigned char touche,int x,int y)
{
  switch (touche)
    {
    case 'p': /* affichage du carre plein */
      glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
      glutPostRedisplay();
      break;
    case 'f': /* affichage en mode fil de fer */
      glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      glutPostRedisplay();
      break;
    case 's' : /* Affichage en mode sommets seuls */
      glPolygonMode(GL_FRONT_AND_BACK,GL_POINT);
      glutPostRedisplay();
      break;
    case 'd':
      glEnable(GL_DEPTH_TEST);
      glutPostRedisplay();
      break;
    case 'D':
      glDisable(GL_DEPTH_TEST);
      glutPostRedisplay();
      break;
    case 'a'://Les faces à l'envers s'affichent en fil de fer
      glPolygonMode(GL_FRONT,GL_FILL);
      glPolygonMode(GL_FRONT,GL_LINE);
      glutPostRedisplay();
    break;
    case 'r' :
        AngleTete+=10.0;
        glutPostRedisplay();
        break;
    case 'R':
        AngleTete+= -10.0;
        glutPostRedisplay();
        break;
    case 'm' :
        AnglePieds+=5.0;
        glutPostRedisplay();
        break;
    case 'M':
        AnglePieds+= -5.0;
        glutPostRedisplay();
        break;
    case 'e':
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glutPostRedisplay();
    break;
    case 'E':
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    glutPostRedisplay();
    break;
    case 'Z':
      ZoomZ = ZoomZ + 0.5;
      glutPostRedisplay();
      break;
    case 'z':
      ZoomZ = ZoomZ - 0.5;
      glutPostRedisplay();
      break;
    case 'c':
		lum = 1;
    break;
    case 'v':
        lum = 2;
    break;
    case 'b':
        lum = 3;
    break;
    case 'n':
        lum = 0;
    break;
    case 'k':
        mat = 1;
		glutPostRedisplay();
        break;
    case 'l':
        mat = 2;
        glutPostRedisplay();
        break;
    case 'L':
        mat = 0;
        glutPostRedisplay();
        break;
    case 'q' : /*la touche 'q' permet de quitter le programme */
      exit(0);
    }
}

void specialKeys(int touche, int x, int y) {
    switch (touche) {
        case GLUT_KEY_UP:
            RotaVertical = RotaVertical +2;
            glutPostRedisplay();
            break;
        case GLUT_KEY_DOWN:
            RotaVertical = RotaVertical -2;
            glutPostRedisplay();
            break;
        case GLUT_KEY_LEFT:
            RotaHorizontale = RotaHorizontale +2;
            glutPostRedisplay();
            break;
        case GLUT_KEY_RIGHT:
            RotaHorizontale = RotaHorizontale -2;
            glutPostRedisplay();
            break;
    }
}

void reshape(int x,int y)
{
  if (x<y)
    glViewport(0,(y-x)/2,x,x);
  else
    glViewport((x-y)/2,0,y,y);
}

void mouse(int button, int state,int x,int y)
{
  /* si on appuie sur le bouton gauche */
  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
  {
    presse = 1; /* le booleen presse passe a 1 (vrai) */
    xold = x; /* on sauvegarde la position de la souris */
    yold=y;
  }
  /* si on relache le bouton gauche */
  if (button == GLUT_LEFT_BUTTON && state == GLUT_UP)
    presse=0; /* le booleen presse passe a 0 (faux) */
}

void mousemotion(int x,int y)
  {
    if (presse) /* si le bouton gauche est presse */
    {
      /* on modifie les angles de rotation de l'objet
	 en fonction de la position actuelle de la souris et de la derniere
	 position sauvegardee */
      anglex=anglex+(x-xold);
      angley=angley+(y-yold);
      glutPostRedisplay(); /* on demande un rafraichissement de l'affichage */
    }

    xold=x; /* sauvegarde des valeurs courante de le position de la souris */
    yold=y;
  }

//méthode pour l'animation automatique des oreilles
void animationOreille()
{
    oreilAngle += oreilRapidite;
    if(oreilAngle > 40 || oreilAngle < 0)
    {
        oreilRapidite = -oreilRapidite;
    }
    glutPostRedisplay();
}

void animationBras()
{
    static float angleRotation = 0.0;
    static float direction = 1.0;
    const float maxAngle = 10.0;  // Angle maximal de rotation

    angleRotation += BrasRapid * direction;

    if (angleRotation >= maxAngle || angleRotation <= -maxAngle) {
        direction *= -1.0;  // Inverse la direction de rotation lorsque l'angle maximal est atteint
    }

    AngleBras = angleRotation;
    glutPostRedisplay();
}

//méthode pour l'animation avec les touches du clavier
void animationTete()
{
    AngleTete += TeteRapid;

    glutPostRedisplay();
}

//méthode pour l'animation avec les touches du clavier
void animationPieds()
{
    AnglePieds += PiedsRapid;

    glutPostRedisplay();
}

//méthode pour charger l'image jouant le rôle de texture
void loadJpegImage(char *fichier)
{
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  FILE *file;
  unsigned char *ligne;

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
#ifdef __WIN32
  if (fopen_s(&file,fichier,"rb") != 0)
    {
      fprintf(stderr,"Erreur : impossible d'ouvrir le fichier texture.jpg\n");
      exit(1);
    }
#elif __GNUC__
  if ((file = fopen(fichier,"rb")) == 0)
    {
      fprintf(stderr,"Erreur : impossible d'ouvrir le fichier texture.jpg\n");
      exit(1);
    }
#endif
  jpeg_stdio_src(&cinfo, file);
  jpeg_read_header(&cinfo, TRUE);

  if ((cinfo.image_width!=256)||(cinfo.image_height!=256)) {
    fprintf(stdout,"Erreur : l'image doit etre de taille 256x256\n");
    exit(1);
  }
  if (cinfo.jpeg_color_space==JCS_GRAYSCALE) {
    fprintf(stdout,"Erreur : l'image doit etre de type RGB\n");
    exit(1);
  }

  jpeg_start_decompress(&cinfo);
  ligne=image;
  while (cinfo.output_scanline<cinfo.output_height)
    {
      ligne=image+3*256*cinfo.output_scanline;
      jpeg_read_scanlines(&cinfo,&ligne,1);
    }
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
}
