from pathlib import Path
import gzip
import zipfile


class FileTool:
    @staticmethod
    def get_project_dir(sub=None, *, relative=False):
        """
        Root directory of PyGWater (the folder that contains the src_* packages and 'data').
        Resolved from the location of this file, so it no longer depends on the current working directory.
        """
        dir_of_project = Path(__file__).resolve().parents[1]

        if relative:
            import os
            result = Path(os.path.relpath(dir_of_project, Path().absolute()))
        else:
            result = dir_of_project

        if sub is not None:
            result /= sub

        return result

    @staticmethod
    def get_data_dir(sub=None):
        """
        External data root, kept outside the code (cf. PyGLDA_v2_external_data): the environment variable
        PYGWATER_DATA if set, otherwise <project>/data. It holds the inputs (LLN/, land_mask/) and the files written
        at run time (res/ for results, temp/ for scratch), see get_res_dir / get_temp_dir.
        """
        import os
        root = os.environ.get('PYGWATER_DATA')
        result = Path(root) if root else FileTool.get_project_dir('data')
        if sub is not None:
            result /= sub
        return result

    @staticmethod
    def get_res_dir(sub=None):
        """results written by PyGWater (templates, masks, shapefiles, figures data): <data root>/res"""
        return FileTool.get_data_dir('res') if sub is None else FileTool.get_data_dir('res') / sub

    @staticmethod
    def get_temp_dir(sub=None):
        """scratch files: <data root>/temp"""
        return FileTool.get_data_dir('temp') if sub is None else FileTool.get_data_dir('temp') / sub

    @staticmethod
    def get_files_in_dir(fp, sub=False):
        assert fp.is_dir, f"{fp} is not a directory."

        file_list = []

        iterlist = list(fp.iterdir())
        for i in range(len(iterlist)):
            if iterlist[i].is_file():
                file_list.append(iterlist[i])

            elif sub:
                file_list += FileTool.get_files_in_dir(iterlist[i], sub=sub)

        return file_list

    @staticmethod
    def un_gz(gz_file_path: Path, target_path: Path = None):

        if target_path is None:
            target_path = gz_file_path.parent / gz_file_path.name.replace('.gz', '')

        g_file = gzip.GzipFile(gz_file_path)
        open(target_path, "wb+").write(g_file.read())
        g_file.close()

    @staticmethod
    def un_zip(zip_file_path: Path, target_path: Path = None):
        if target_path is None:
            target_path = zip_file_path.parent / zip_file_path.name.replace('.zip', '')

        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            zip_ref.extractall(target_path)
