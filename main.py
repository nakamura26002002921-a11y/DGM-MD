
pip install huggingface_hub

python -c "from huggingface_hub import hf_hub_download; hf_hub_download(repo_id='bartowski/gemma-2-2b-it-GGUF', filename='gemma-2-2b-it-Q4_K_M.gguf', local_dir='/home/models', local_dir_use_symlinks=False)"
