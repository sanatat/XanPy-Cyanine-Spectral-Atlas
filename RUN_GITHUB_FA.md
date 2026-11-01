# راهنمای فارسی آماده‌سازی و انتشار در GitHub

این مخزن داده‌های خام مقاله را **دوباره منتشر نمی‌کند**؛ برای بازتولید، چهار فایل ZIP اصلی باید در `data/raw/` قرار بگیرند، یا با دستور زیر از Zenodo دریافت شوند.

## مرحله اول: از ZIP به یک پوشه پروژه

```bash
cd ~/Downloads
unzip -o XanPy-Cyanine-Atlas-GitHub.zip
cd XanPy-Cyanine-Atlas-GitHub
```

## مرحله دوم: نصب و بررسی (استفاده از محیط مجازی قبلی)

```bash
~/Downloads/.xanpy-venv/bin/python -m pip install -e '.[test]'
~/Downloads/.xanpy-venv/bin/python -m pytest -q
```

## مرحله سوم: اجرای تحلیل روی فایل‌هایی که قبلاً دانلود کرده‌ای

```bash
~/Downloads/.xanpy-venv/bin/python -m xanpy_cyanine.analysis \
  --zip-dir ~/Downloads/XanPy-Stage19-Independent-Validation/raw_zips \
  --out-dir results
```

یا اگر فایل‌های قبلی را نداری و می‌خواهی مستقیماً از Zenodo دانلود کنی:

```bash
~/Downloads/.xanpy-venv/bin/python scripts/download_zenodo.py
~/Downloads/.xanpy-venv/bin/python -m xanpy_cyanine.analysis --zip-dir data/raw --out-dir results
```

## مرحله چهارم: ساخت مخزن جدید

1. در سایت https://github.com/new یک مخزن **خالی** با نام `XanPy-Cyanine-Spectral-Atlas` بساز؛ گزینه ساخت README جدید را فعال نکن.
2. در ترمینال داخل پوشه همین پروژه:

```bash
git init
git add .
git status --short
git commit -m "Initial release: reproducible cyanine spectral-shape reanalysis"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/XanPy-Cyanine-Spectral-Atlas.git
git push -u origin main
```

به‌جای `YOUR_USERNAME` نام کاربری GitHub خودت را بگذار. برای `git push` به احراز هویت GitHub نیاز داری؛ ممکن است از SSH یا GitHub CLI هم استفاده کنی.

**قبل از push بررسی کن:** فایل‌های `Cy5-7.zip` و بقیه فایل‌های خام نباید در خروجی `git status --short` دیده شوند. داده‌های حجیم/خام داخل `.gitignore` هستند.

## چه نتایجی را در README برجسته کرده‌ایم؟

- تغییر جذب Cy5 حدود **−۰٫۳۸ nm** است؛ با وجود این، پهنای ظاهری پوش طیف جذب حدود **۳۹٫۲٪ کمتر** می‌شود.
- تغییر جذب Cy7 حدود **−۴٫۶۱ nm** و تغییر نشر حدود **−۱۷٫۱۱ nm** است.
- این یافته‌ها *بازتحلیل داده‌های آرشیوی پژوهشگران اصلی* هستند. بدون اسکن تکراری نمی‌توانیم ادعای معناداری آماری یا تعیین سازوکار کنیم.

گزارش کامل انگلیسی: `docs/RESULTS.md` و `docs/METHODS.md`.
