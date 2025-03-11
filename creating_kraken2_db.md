# Instructions for recreating/modifying ADAR's Kraken2 database

1. **Ensure you have at least 50GB+ of storage available**.
2. Download and install [Kraken2](https://github.com/DerrickWood/kraken2)
3. Create a new directory for subsequent steps, and `cd` into it.
4. Run `kraken2-build --download-taxonomy --db $DB_NAME`; Bear in mind that your choice for `$DB_NAME` will be used to create a directory.
5. If you want to start with the sequences from our own reference database, download `assets/ref.fa` from this repository. You'll need to create individual fastas from each entry, which you can achieve in python with something like this: 

    ```python
    import os

    if not os.path.exists("refs"):
        os.mkdir("refs")

    with open("ref.fa", "r") as inf:
        lines = inf.readlines()

    for i in range(0, len(lines), 2):
        header = lines[i]
        name = header[1:]
        seq = lines[i + 1]
        outfile = os.path.join("refs", f"{name.split(" ")[0]}.fa")

        with open(outfile, "w") as outf:
            outf.write(header)
            outf.write(seq)
    ```

    Each fasta will be found inside `./refs/`.

6. download any other fastas you may need, and store them in the same folder as the ref.fa references, if you want to use them.
7. add each entry to your fresh kraken2 db. This can be done at once for all files within a directory (`refs` in this example):
```bash
for file in refs/*.fa; do
	kraken2-build --add-to-library $file --db $DB_NAME
done
```
8. Once you've added all entries to the database, build it with this command:
```bash
kraken2-build --build -db $DB_NAME
```
9. The finished database will be a directory called `$DB_NAME`. Now you can run adar with `--kraken2_db $PATH_TO_DB_NAME`.

> [!NOTE]
> If you don't want to add any more references to your Kraken2 database and you want to save ~44GB of space, delete the `taxonomy` directory within the database. If you ever need to add more sequences to the database, you'll have to redownload it.

