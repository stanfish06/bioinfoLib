# Copyright 2025 Zhiyuan Yu (Heemskerk's lab, University of Michigan)

import mysql.connector
from fuzzywuzzy import fuzz


def query_gene(
    g: str, db_name: str, db_type: str = "core", tb_name: str = "gene", port: int = 3306
):
    cnx = mysql.connector.connect(
        user="anonymous", host="ensembldb.ensembl.org", port=port
    )
    print(f"server version: {cnx.server_info}")
    c = cnx.cursor()
    c.execute("SHOW DATABASES")
    databases = c.fetchall()
    databases = [
        db for db in databases if fuzz.ratio(db[0], db_name) > 50 and (db_type in db[0])
    ]
    print("selected databases:")
    databases = list(map(lambda x: (print(x[0]), x[0])[1], databases))
    c.execute(f"USE {databases[0]}")
    c.execute(f"DESCRIBE {tb_name}")
    out = c.fetchall()
    print(out)
    c.execute(f"SELECT * FROM {tb_name} LIMIT 50")
    out = c.fetchall()
    print(out)
    cnx.close()


if __name__ == "__main__":
    query_gene("SOX2", "homo sapien", "core")
