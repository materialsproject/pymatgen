#ifdef STAR_SCANNER
/* The defines give the position of the token in this
array */
char * tokens[] = {"END","LBLOCK","GLOBAL","STOP",
                 "save_heading","save_end",
		 "data_heading","data_name",
		 "start_sc_line","sc_line_of_text",
		 "end_sc_line","data_value_1",
		 "(error)"};
char * input_string;     /* where flex gets input */
size_t string_pos;          /* current position */
size_t in_string_len;       /* total length */
int * token_list;           /* list of tokens */
int * line_no_list;         /* list of token positions */
char ** value_list;        /* list of values */
size_t alloc_mem;           /* How much allocated */
size_t current_len;         /* Length of list */
#define MEM_ALLOC_SIZE 4192 /* Allocation block size */
extern int star_scanner(void);
extern void star_clear(void);
extern char * yytext;
extern size_t yyleng;
extern size_t yylineno;
#else
extern char * input_string;
extern size_t string_pos;
extern size_t in_string_len;
#endif

/* defines to index into the token list */

#define DEND 0
#define DLBLOCK 1
#define DGLOBAL 2
#define DSTOP 3
#define DSAVE_HEADING 4
#define DSAVE_END 5
#define DDATA_HEADING 6
#define DDATA_NAME 7
#define DSTART_SC_LINE 8
#define DSC_LINE_OF_TEXT 9
#define DEND_SC_LINE 10 
#define DDATA_VALUE_1 11
#define DERROR 12
