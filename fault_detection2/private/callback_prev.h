











/*80*************************************************************************************/	
//
// Autres
//
/*80*************************************************************************************/

#define attDataAddSurface3d	\
				RET_PVOID(long id_object, \
									double x1, double y1, double z1, \
									double x2, double y2, double z2, \
									double x3, double y3, double z3, \
									double x4, double y4, double z4, \
									double r, double g, double b, \
									double * texture, long width_texture, long height_texture) \
				ATT_CALL("attDataAddSurface3d")))

#define attDataGetPositionFromScreen \
				RET_VOID(long id_object, long x_screen, long y_screen, double * x, double * y, double * z) \
				ATT_CALL("attDataGetPositionFromScreen")))

#define attDataAnimateBegin \
				RET_VOID(long id_object) \
				ATT_CALL("attDataAnimateBegin")))

#define attDataAnimateEnd \
				RET_VOID(long id_object) \
				ATT_CALL("attDataAnimateEnd")))

#define attDataCallbackGetMainHwnd \
				RET_PVOID() \
				ATT_CALL("attDataCallbackGetMainHwnd")))

#define attDataCallbackGetHwnd \
				RET_PVOID(long id_object) \
				ATT_CALL("attDataCallbackGetHwnd")))

#define attDataGetSpace \
				RET_LONG(long id_object, double * offset, double * scale) \
				ATT_CALL("attDataGetSpace")))

#define attDataGetView3D \
				RET_LONG(long id_object) \
				ATT_CALL("attDataGetView3D")))

#define attDataAddArrowRGB \
				RET_VOID(long id, double x1, double y1, double x2, double y2, double norm_arrow, const char * name,\
								 double r, double g, double b, long thickness) \
				ATT_CALL("attDataAddArrowRGB")))

#define attFileInit	\
				RET_PVOID(long hwnd, char * title) \
				ATT_CALL("attFileInit")))

#define attFileSelectorOpen	\
				RET_LONG(void * data) \
				ATT_CALL("attFileSelectorOpen")))

#define attFileGetName \
				RET_PCHAR(void * data) \
				ATT_CALL("attFileGetName")))

#define attFileSetName \
				RET_VOID(void * data, const char * name) \
				ATT_CALL("attFileSetName")))

#define attFileDestroy \
				RET_VOID(void * data) \
				ATT_CALL("attFileDestroy")))

#define attExec	\
				RET_VOID(char * function, ...) \
				ATT_CALL("attExec")))


/*80*************************************************************************************/	
//
// Obsolete
//
/*80*************************************************************************************/

typedef struct _ATT_DESC_DIALOG
{
	char title[200];	
	long x, y, dx, dy;
	long style;
	long font_size;
	char font_name[200];
	ndfunc mainproc;	
	struct
	{
		long type;
		ndfunc proc;
	} msg[200];	
	struct		
	{
		long type;
		char name[200];
		long id;
		long x, y, dx, dy;
	} entry[500];
} ATT_DESC_DIALOG;


/* BOITES DE DIALOGUES */
#define DIALOG_NO								0
#define DIALOG_GROUPBOX					1
#define DIALOG_FINGROUPBOX			2
#define DIALOG_LTEXT						3
#define DIALOG_RTEXT						4
#define DIALOG_CTEXT						5
#define DIALOG_EDIT							6
#define DIALOG_RADIO						7
#define	DIALOG_CHECKBOX					8
#define DIALOG_PUSHBUTTON				9
#define DIALOG_OWNER						10
#define DIALOG_OWNER_EX					11
#define DIALOG_TRACKBAR					12
#define END_DIALOG						DIALOG_NO, "", 0, 0, 0, 0, 0
#define DIALOG_END_GROUPBOX		DIALOG_FINGROUPBOX, "", 0, 0, 0, 0, 0,
#define MSG_NO						0
#define MSG_INIT					1
#define MSG_CMD						2
#define MSG_CLOSE					3
#define MSG_INITWITHDATA	4
#define END_MSG						MSG_NO, (attfunc)NULL
#define NO_MSG						{ END_MSG },
#define FONT_8						8, "MS Sans Serif",
#define STYLE_MIN					0
#define STYLE_NORMAL			1
#define STYLE_FIN					2
#define STYLE_3D					3

#define attDialogRunFromDESC \
				RET_VOID(ATT_DESC_DIALOG * dialog) \
				ATT_CALL("attDialogRunFromDESC")))

#define attDialogWithDataRunFromDESC \
				RET_VOID(ATT_DESC_DIALOG * dialog, void * data) \
				ATT_CALL("attDialogWithDataRunFromDESC")))

#define attDialogModalRunFromDESC	\
				RET_LONG(ATT_DESC_DIALOG * dialog) \
				ATT_CALL("attDialogModalRunFromDESC")))

#define attDialogModalWithDataRunFromDESC	\
				RET_LONG(ATT_DESC_DIALOG * dialog, void * data) \
				ATT_CALL("attDialogModalWithDataRunFromDESC")))

#define attDialogEnd \
				RET_VOID(long hwnd, long end) \
				ATT_CALL("attDialogEnd")))

#define attWindowDestroy \
				RET_LONG(long hwnd) \
				ATT_CALL("attWindowDestroy")))

#define attWindowEnableItem	\
				RET_VOID(long hwnd, long id, long state) \
				ATT_CALL("attWindowEnableItem")))

#define attWindowCheckRadioButton	\
				RET_VOID(long hwnd, long id1, long id2, long id_active) \
				ATT_CALL("attWindowCheckRadioButton")))

#define attWindowCheckButton \
				RET_VOID(long hwnd, long id, long state) \
				ATT_CALL("attWindowCheckButton")))

#define attWindowModifyItemReal	\
				RET_VOID(long hwnd, long id, double * variable, long sign) \
				ATT_CALL("attWindowModifyItemReal")))

#define attWindowModifyItemRealMinMax	\
				RET_VOID(long hwnd, long id, double * variable, double min, double max) \
				ATT_CALL("attWindowModifyItemRealMinMax")))

#define attWindowModifyItemInt \
				RET_VOID(long hwnd, long id, long * variable, long even, long sign) \
				ATT_CALL("attWindowModifyItemInt")))

#define attWindowModifyItemIntMinMax \
				RET_VOID(long hwnd, long id, long * variable, long min, long max) \
				ATT_CALL("attWindowModifyItemIntMinMax")))

#define attWindowSetItemReal \
				RET_VOID(long hwnd, long id, double l) \
				ATT_CALL("attWindowSetItemReal")))

#define attWindowSetItemReal2	\
				RET_VOID(long hwnd, long id, double l) \
				ATT_CALL("attWindowSetItemReal2")))

#define attWindowSetItemInt	\
				RET_VOID(long hwnd, long id, long variable) \
				ATT_CALL("attWindowSetItemInt")))

#define attWindowGetItemText \
				RET_VOID(long hwnd, long id, char * string, long size) \
				ATT_CALL("attWindowGetItemText")))

#define attWindowSetItemText \
				RET_VOID(long hwnd, long id, char * string) \
				ATT_CALL("attWindowSetItemText")))

/* Types */
typedef int BOOL;		

/* Identificateurs */
#if !defined (TRUE)
	#define TRUE	1
#endif
#if !defined (FALSE)
	#define FALSE	0
#endif
#if !defined (NULL)
	#define NULL 0
#endif
#if !defined (PI_DOUBLE)
	#define PI_DOUBLE  6.28318530718
#endif
#if !defined (PI)
	#define PI         3.14159265359	
#endif
#if !defined (DEMI_PI)
	#define DEMI_PI		 1.57079632679	
#endif

/* Fonctions Mathématiques */
#if !defined (MIN)
	#define MIN(x,y)		( ( x >= y ) ? y : x )
#endif
#if !defined (MAX)
	#define MAX(x,y)		( ( x >= y ) ? x : y )
#endif
#if !defined (SGN)
	#define SGN(x)			( ( x >= 0 ) ? (1) : (-1) )
#endif	
#if !defined (RDOUBLE)
	#define RDOUBLE(x)	( floor(x+0.5) )
#endif
#if !defined (RLONG)
	#define RLONG(x)		( (long)floor(x+0.5) )
#endif
#if !defined (CLONG)
	#define CLONG(x)		( (long)ceil(x+0.5) )
#endif

#define DEGFROMRAD	 57.29577951307
#define RADTODEG(x)	( x*DEGFROMRAD)
#define DEGTORAD(x)	( x/DEGFROMRAD)

#define LIST_INC(x)						x = (x)->next;
#define LIST_DEC(x)						x = (x)->prev;
#define LIST_ADDR_INC(x)			x = &((*x)->next);
#define LIST_LAST(list)				while ( list->next != NULL ) LIST_INC(list)
#define IS_LIST_NEXT(list)		( (list)->next != NULL )
#define IS_LIST_PREV(list)		( (list)->prev != NULL )
#define LIST_NEXT(x)					( (x)->next )
#define LIST_PREV(x)					( (x)->prev )
#define LIST_DETACH(x)				(x)->prev = (x)->next = NULL;

#define LIST_ATTACH(list_element_prev, list_attach) \
	{ \
		if ( list_element_prev != NULL ) (list_element_prev)->next = list_attach; \
		if ( list_attach != NULL ) (list_attach)->precedent = list_element_prev; \
	}
#define LIST_INSERT(list_element_prev, list_element_insert) \
	{ \
		(list_element_insert)->prev = list_element_prev; \
		(list_element_insert)->next = LIST_NEXT(list_element_prev); \
		(list_element_prev)->next = list_element_insert; \
		if ( LIST_NEXT(list_element_prev) ) (list_element_prev)->next->prev = list_element_insert; \
	}
#define LIST_INSERT_BEGIN(list_begin, list_element) \
	{ \
		(list_element)->next = list_begin; \
		if ( list_begin != NULL ) (list_begin)->prev = list_element; \
		list_begin = list_element; \
	}
#define LIST_INSERT_AFTER_BEGIN(list_begin, list_element) \
	{ \
		if ( list_begin == NULL ) list_begin = list_element; \
		else \
			{ \
				(list_element)->prev = list_begin; \
				(list_element)->next = LIST_NEXT(list_begin); \
				if ( IS_LIST_PREV(list_element) ) (list_element)->prev->next = list_element; \
				if ( IS_LIST_NEXT(list_element) ) (list_element)->next->prev = list_element; \
			} \
	}
#define LIST_NEW(list_begin, list_element, list_size) \
	{ \
		list_element = (list_size *) calloc (1, sizeof(list_size)); \
		LIST_INSERT_BEGIN(list_begin, list_element) \
	}
#define LIST_CUT(list_element) \
	{ \
		if ( IS_LIST_NEXT(list_element) ) (list_element)->next->prev = LIST_PREV(list_element); \
		if ( IS_LIST_PREV(list_element) ) (list_element)->prev->next = LIST_NEXT(list_element); \
	}
#define LIST_DELETE(list_begin, list_element) \
	{ \
		LIST_CUT(list_element) \
		if ( list_begin == list_element ) list_begin = LIST_NEXT(list_element); \
		free (list_element); \
	}
