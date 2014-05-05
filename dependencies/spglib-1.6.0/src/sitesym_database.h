/* sitesym_database.h */
/* Copyright (C) 2011 Atsushi Togo */

#ifndef __sitesym_database_H__
#define __sitesym_database_H__

int ssmdb_get_coordinate( int rot[3][3],
			  double trans[3],
			  const int index );
void ssmdb_get_wyckoff_indices( int indices[2], const int index );

#endif
