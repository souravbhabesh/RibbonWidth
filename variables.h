#define PERIOD 10000
#define NMAX 20000
#define a 1.0
#define NXMAX 201
#define MAXFRAMES 20001
#define MAXRUN 201
#define MAXPARTICLETYPE 10

extern int N,Nb,Nd,bondGroup[NMAX*2],dihedralGroup[NMAX*4];
//N:#particles, Nb:#bonds, Nd:#dihedrals
extern float position[NMAX*3];
extern uint32_t particleID[NMAX];
extern char particleType[MAXPARTICLETYPE][2];

extern int NX,NY,RUN,STEPS,LEN,FRAMES,JK_BIN_COUNT;
extern double KAPPA,EPSILON;

//Jack Knife blocking variables
extern double topx[NXMAX][MAXRUN][MAXFRAMES];
extern double topy[NXMAX][MAXRUN][MAXFRAMES];
extern double bottomx[NXMAX][MAXRUN][MAXFRAMES];
extern double bottomy[NXMAX][MAXRUN][MAXFRAMES];

// Ribbon width 
extern double framewidth[NXMAX][MAXRUN][MAXFRAMES];
extern double rwidth[NXMAX][MAXRUN];

