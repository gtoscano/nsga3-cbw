/* Routines for storing population data into files */

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <arrow/type_fwd.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <parquet/exception.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>

#include "global.h"
#include "rand.h"


/* Function to print the information of a population in a file */
void report_pop(population *pop, FILE *fpt)
{
    int i, j, k;
    for (i = 0; i < popsize; i++)
    {
        for (j = 0; j < nobj; j++)
        {
            fprintf(fpt, "%e\t", pop->ind[i].obj[j]);
        }
        if (ncon != 0)
        {
            /*
            for (j = 0; j < ncon; j++)
            {
                fprintf(fpt, "%e\t", pop->ind[i].constr[j]);
            }
            */
        }
        if (nreal != 0)
        {
            /*
            for (j = 0; j < nreal; j++)
            {
                fprintf(fpt, "%e\t", pop->ind[i].xreal[j]);
            }
            */
        }
        if (nbin != 0)
        {
            for (j = 0; j < nbin; j++)
            {
                for (k = 0; k < nbits[j]; k++)
                {
                    fprintf(fpt, "%d\t", pop->ind[i].gene[j][k]);
                }
            }
        }
        fprintf(fpt, "%e\t", pop->ind[i].constr_violation);
        fprintf(fpt, "%d\n", pop->ind[i].rank);
    }
    return;
}

/* Function to print the information of feasible and non-dominated population in a file */

// #1 Write out the data as a Parquet file
void write_parquet_file(const arrow::Table& table) {
  std::shared_ptr<arrow::io::FileOutputStream> outfile;
  PARQUET_ASSIGN_OR_THROW(
      outfile, arrow::io::FileOutputStream::Open("parquet-arrow-example.parquet"));
  // The last argument to the function call is the size of the RowGroup in
  // the parquet file. Normally you would choose this to be rather large but
  // for the example, we use a small value to have multiple RowGroups.
  PARQUET_THROW_NOT_OK(
      parquet::arrow::WriteTable(table, arrow::default_memory_pool(), outfile, (int) (nobj+nreal)/4));
}


void read_single_column() {
  std::shared_ptr<arrow::io::ReadableFile> infile;
  PARQUET_ASSIGN_OR_THROW(infile,
                          arrow::io::ReadableFile::Open("parquet-arrow-example.parquet",
                                                        arrow::default_memory_pool()));

  std::unique_ptr<parquet::arrow::FileReader> reader;
  PARQUET_THROW_NOT_OK(
      parquet::arrow::OpenFile(infile, arrow::default_memory_pool(), &reader));
  std::shared_ptr<arrow::ChunkedArray> array;
  PARQUET_THROW_NOT_OK(reader->ReadColumn(4, &array));
  PARQUET_THROW_NOT_OK(arrow::PrettyPrint(*array, 4, &std::cout));
}

void report_feasible(population *pop, FILE *fpt)
{
    int i, j, k;
    int ncols = nobj + nreal;
    int nrows = 0;
    for (i = 0; i < popsize; i++)
    {
        if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank == 1)
        {
          ++nrows;
        }
    }

    std::vector< std::vector<double> > mat(ncols, std::vector<double>(nrows));
    std::vector< double> tmp_row (ncols);
    int pop_counter = 0;

    std::vector<std::shared_ptr<arrow::Field> > schema_labels_vec;
    for (i = 0; i < popsize; i++)
    {
        if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank == 1)
        {
            std::copy(pop->ind[i].obj, pop->ind[i].obj+nobj, tmp_row.begin());

            for (j = 0; j < nobj; j++)
            {
                fprintf(fpt, "%e\t", pop->ind[i].obj[j]);
                mat[j][pop_counter] = pop->ind[i].obj[j];
                if (pop_counter == 0) {

                  char tmp_str[20];
                  sprintf(tmp_str, "obj_%i",j);
                  schema_labels_vec.push_back(arrow::field(tmp_str, arrow::float64()));
                }

            }
            /*
            if (ncon != 0)
            {
                for (j = 0; j < ncon; j++)
                {
                    fprintf(fpt, "%e\t", pop->ind[i].constr[j]);
                }
            }
            */
            if (nreal != 0)
            {
                std::copy(pop->ind[i].xreal, pop->ind[i].xreal+nreal, tmp_row.begin() + nobj);
                for (j = 0; j < nreal; j++)
                {
                    fprintf(fpt, "%e\t", pop->ind[i].xreal[j]);
                    mat[nobj+j][pop_counter] = pop->ind[i].xreal[j];
                    if (pop_counter == 0) {
                      char tmp_str[20];
                      sprintf(tmp_str, "x_%i",j);

                      schema_labels_vec.push_back(arrow::field(tmp_str, arrow::float64()));
                    }
                }
            }
            if (nbin != 0)
            {
                for (j = 0; j < nbin; j++)
                {
                    for (k = 0; k < nbits[j]; k++)
                    {
                        fprintf(fpt, "%d\t", pop->ind[i].gene[j][k]);
                    }
                }
            }
            ++pop_counter;
            //fprintf(fpt, "%e\t", pop->ind[i].constr_violation);
            //fprintf(fpt, "%d\n", pop->ind[i].rank);
            fprintf(fpt,"\n");

        }

    }
    return;

    std::shared_ptr<arrow::Schema> schema = arrow::schema(schema_labels_vec);
    std::vector< std::shared_ptr<arrow::Array> > schema_values_vec;
    for (auto col : mat){
      /*
      arrow::DoubleBuilder double_builder;
      PARQUET_THROW_NOT_OK(double_builder.AppendValues(col));
      std::shared_ptr<arrow::Array> double_array;
      PARQUET_THROW_NOT_OK(double_builder.Finish(&double_array));
      schema_values_vec.push_back(double_array);
      */

      arrow::DoubleBuilder double_builder;
      PARQUET_THROW_NOT_OK(double_builder.AppendValues(col));
      std::shared_ptr<arrow::Array> double_array;
      PARQUET_THROW_NOT_OK(double_builder.Finish(&double_array));
      schema_values_vec.push_back(double_array);
    }

  std::shared_ptr<arrow::Table> table = arrow::Table::Make(schema, schema_values_vec);
  //write_parquet_file(*table);
    //try {
      std::shared_ptr<arrow::io::FileOutputStream> outfile;
      PARQUET_ASSIGN_OR_THROW(
        outfile,
        arrow::io::FileOutputStream::Open("pareto_set.parquet"));

      PARQUET_THROW_NOT_OK(
        parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, (int)ncols/4));

      //PARQUET_ASSIGN_OR_THROW(
      //    outfile, arrow::io::FileOutputStream::Open("parquet-arrow-example.parquet"));
      // The last argument to the function call is the size of the RowGroup in
      // the parquet file. Normally you would choose this to be rather large but
      // for the example, we use a small value to have multiple RowGroups.
      //PARQUET_THROW_NOT_OK(
      //  parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, (int)ncols/4));
    /*} catch (exception &ex) {
      std::cout<<ex.what()<<std::endl;
    }*/
    std::cout<<std::endl;
    read_single_column();
    return;
}
