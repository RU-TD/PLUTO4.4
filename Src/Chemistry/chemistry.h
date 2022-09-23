/* ############################################################

     FILE:     chemistry.h

     PURPOSE:  contains shared definitions with scope
               limited to the chemistry module ONLY

   ############################################################ */

/* ##############################################################

                   P R O T O T Y P I N G

   ############################################################## */

void   initialize_Microphysics(Grid *grid);
void   cleanup_Microphysics();
void   Chemistry(Data_Arr v, double dt, Grid *grid);
void   find_CommunicationNeighbour(int current_rank, LocalDomainInfo *domain_info_array, int nproc, CommunicationNeighbour* cn);
void   find_CommunicationNeighbours(Grid *grid);
void   calculate_Attenuation_perDomain(Data_Arr v, Grid *grid);
void   calculate_ColumnDensity_perDomain(Data_Arr v, Grid *grid, int val);
void   calculate_ColumnDensity(Data_Arr v, Grid *grid);
void   calculate_Attenuation(Data_Arr v, Grid *grid);
