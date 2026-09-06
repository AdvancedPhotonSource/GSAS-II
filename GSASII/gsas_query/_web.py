"""Entry point for `gsas2-query-web` console script."""

import os


def main():
    import uvicorn
    from .app import app

    host = os.environ.get("HOST", "0.0.0.0")
    port = int(os.environ.get("PORT", "8000"))
    uvicorn.run(app, host=host, port=port)


if __name__ == "__main__":
    main()
