
pip install huggingface_hub

python -c "from huggingface_hub import hf_hub_download; hf_hub_download(repo_id='unsloth/gemma-4-E4B-it-GGUF', filename='gemma-4-E4B-it-Q4_K_M.gguf', local_dir='/home')"
