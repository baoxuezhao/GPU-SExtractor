/*
 globals.h

 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *
 *	Part of:	SExtractor
 *
 *	Author:		E.BERTIN (IAP)
 *
 *	Contents:	global declarations.
 *
 *	Last modify:	11/05/2008
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */

#include	"types.h"

/*----------------------- miscellaneous variables ---------------------------*/

extern sexcatstruct thecat;
extern picstruct 	thefield1, thefield2, thewfield1, thewfield2;
extern objstruct 	flagobj;
extern obj2struct 	flagobj2;
extern objstruct	outobj;
extern obj2struct 	outobj2;
extern float ctg[37], stg[37];
extern char gstr[MAXCHAR];

extern prefstruct prefs;
extern brainstruct	*brain;

extern keystruct	objkey[];
/*------------------------------- functions ---------------------------------*/
extern void alloccatparams(void);
extern void allocparcelout(void);
extern void analyse(picstruct *, picstruct *, int, objliststruct *);
extern void blankit(char *, int);
extern void endcat(char *error);
extern void reendcat(void);
extern void changecatparamarrays(char *keyword, int *axisn, int naxis);
extern void closecheck(void);
extern void copydata(picstruct *, int, int);
extern void dumpparams(void);
extern void endfield(picstruct *);
extern void endobject(picstruct *,picstruct *, picstruct *, picstruct *, objstruct *);
extern void examineiso(picstruct *, picstruct *, objstruct *, pliststruct *);
extern void flagcleancrowded(int, objliststruct *);
extern void freeparcelout(void);
extern void getnnw(void), initcat(void);
extern void reinitcat(picstruct *);
extern void initglob(void);
extern void makeit(void);
extern void mergeobject(objstruct *, objstruct *);
extern void neurinit(void);
extern void neurclose(void);
extern void neurresp(double *, double *);
extern void preanalyse(int, objliststruct *, int);
extern void readcatparams(char *);
extern void readdata(picstruct *, PIXTYPE *, int);
extern void readidata(picstruct *, FLAGTYPE *, int);
extern void readimagehead(picstruct *);
extern void readprefs(char *, char **, char **, int);
extern void scanimage(picstruct *, picstruct **, int);
extern void sexcircle(PIXTYPE *bmp, int, int, double, double, double, PIXTYPE);
extern void sexdraw(PIXTYPE *bmp, int, int, double, double, PIXTYPE);
extern void sexellips(PIXTYPE *bmp, int, int,
				double, double, double, double, double, PIXTYPE, int);
extern void sexmove(double, double);
extern void updateparamflags(void);
extern void useprefs(void);
extern void writecat(objstruct *);
extern void write_error(char *msg1, char *msg2);
extern void write_vo_fields(FILE *file);

extern float hmedian(float *, int);

extern int addobj(int, objliststruct *, objliststruct *);
extern int belong(int, objliststruct *, int, objliststruct *);
extern int gatherup(objliststruct *, objliststruct *);
extern int parcelout(objliststruct *, objliststruct *);

extern void *loadstrip(picstruct *);

extern char *readfitshead(FILE *, char *, int *);

extern picstruct *inheritfield(picstruct *infield, int flags);
extern picstruct *newfield(char *, int, int);

