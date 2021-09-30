#ifndef MLW_H_IRRADIATION
#define MLW_H_IRRADIATION

/* structs */
typedef struct LOCAL_DOMAIN_INFO
{
        int x1_begin;
        int x1_end;

        int x2_begin;
        int x2_end;

        int x3_begin;
        int x3_end;

} LocalDomainInfo;

typedef struct COMMUNICATION_NEIGHBOUR
{
        int receive_rank;
        int send_rank;
} CommunicationNeighbour;

typedef struct IRRADIATION_DATA
{
        CommunicationNeighbour neighbour;
        double***       S;
        double***       Tdust;
        double***       Tgas;
        double* data_buffer;
        double* column_density_offset;
} IrradiationData;

extern IrradiationData irradiation;


/* prototypes */
void GetDomainInfo(Grid *grid);
void CDInit(Grid *grid);
void CDFinalise(void);

void FindCommunicationNeighbours(Grid *grid);
void FindCommunicationNeighbour(int current_rank, LocalDomainInfo *domain_info_array, int nproc, CommunicationNeighbour* cn);
void CalculateS(Grid *grid, const Data* data);
void UpdateGasTemp(Grid *grid, const Data* data);
double CalculateT(double column, double rho, double Tdust, double dist);
double read_table_temperature(double column, double xi);

#endif