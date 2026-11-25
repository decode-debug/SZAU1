from pathlib import Path
import subprocess
import tempfile


def main():
    projekt = Path(__file__).resolve().parent
    zrodla = projekt / "pust_szablon_sprawozdania_utf8"
    raport = projekt / "sprawozdanie_Wrobel_Zielinski.pdf"
    with tempfile.TemporaryDirectory(dir=projekt) as katalog:
        polecenie = [
            "pdflatex",
            "-interaction=nonstopmode",
            "-halt-on-error",
            f"-output-directory={katalog}",
            f"-jobname={raport.stem}",
            "sprawozdanie_szablon.tex",
        ]
        for _ in range(2):
            subprocess.run(polecenie, cwd=zrodla, check=True)
        (Path(katalog) / raport.name).replace(raport)
    print(raport)


if __name__ == "__main__":
    main()
