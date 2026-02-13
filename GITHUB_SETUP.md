# Push this folder to a new GitHub repo

This folder is already a git repo with one commit. To put it on GitHub:

1. **Create the repo on GitHub**
   - Go to https://github.com/new
   - Repository name: **Meta_AML_Shiny**
   - Public, no README/license (this folder has its own)
   - Click **Create repository**

2. **Add remote and push** (in a terminal):
   ```bash
   cd /Users/bbenard/Desktop/Meta_AML_Shiny
   git remote add origin https://github.com/YOUR_USERNAME/Meta_AML_Shiny.git
   git branch -M main
   git push -u origin main
   ```
   Replace `YOUR_USERNAME` with your GitHub username.

If you use SSH:
   ```bash
   git remote add origin git@github.com:YOUR_USERNAME/Meta_AML_Shiny.git
   git push -u origin main
   ```

After pushing, you can delete this file or keep it for reference.
