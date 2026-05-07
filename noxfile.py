import nox

_py_versions = range(11, 15)


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test(session):
    session.install("-e", ".", "--group", "test")
    session.chdir("tests")
    session.run(
        "pytest",
        "-s",
        "-x",
        *session.posargs,  # propagates sys.argv to pytest
    )


@nox.session(python=["3.14"])
def testcov(session):
    session.install("-e", ".", "--group", "test")
    session.chdir("tests")
    session.run(
        "pytest",
        "--cov-report",
        "html",
        "--cov",
        "mutation_motif",
    )


@nox.session(python=["3.14"])
def fmt(session):
    session.install("ruff")
    session.run("ruff", "check", "--fix-only", ".")
    session.run("ruff", "format", ".")
