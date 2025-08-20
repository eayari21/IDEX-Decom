import undetected_chromedriver as uc
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time

# === Credentials and Filters ===
USERNAME = "eayari"
PASSWORD = "xukco8-jujpax-fAkxoz"
START_TIME = "2025-093T05:11:55"
END_TIME = "2025-140T05:11:55"
KEYWORD = "ACQ"

# === Launch Browser ===
print("🔧 Launching browser...")
driver = uc.Chrome()
wait = WebDriverWait(driver, 30)

try:
    # 1. Login
    print("🔐 Opening login page...")
    driver.get("https://lasp.colorado.edu/ops/imap/")
    wait.until(EC.presence_of_element_located((By.NAME, "username")))
    print("✅ Login form found.")
    driver.find_element(By.NAME, "username").send_keys(USERNAME)
    driver.find_element(By.NAME, "password").send_keys(PASSWORD + Keys.RETURN)
    print("🔓 Logged in successfully.")

    # 2. Navigate directly to Event Messages
    time.sleep(2)
    print("🌐 Navigating to Event Messages page...")
    driver.get("https://lasp.colorado.edu/ops/imap/flight-ops/event-messages")
    time.sleep(5)
    driver.save_screenshot("after_nav.png")
    print("🖼️ Screenshot saved after navigation.")

    # 3. Find all text inputs
    inputs = driver.find_elements(By.CSS_SELECTOR, "input[type='text']")
    if len(inputs) < 3:
        raise Exception("Expected at least 3 text inputs (start, end, keyword).")
    start_input, end_input, keyword_input = inputs[0], inputs[1], inputs[2]

    # 4. Fill in time range and keyword
    start_input.clear()
    start_input.send_keys(START_TIME)
    end_input.clear()
    end_input.send_keys(END_TIME)
    keyword_input.clear()
    keyword_input.send_keys(KEYWORD)
    print(f"📅 Dates set: {START_TIME} to {END_TIME}")
    print(f"🔑 Keyword set: {KEYWORD}")

    # 5. Click "Include" radio
    include_radio = driver.find_element(By.XPATH, "//label[contains(text(),'Include')]/preceding-sibling::div//input")
    if not include_radio.is_selected():
        driver.execute_script("arguments[0].click();", include_radio)
    print("☑️ 'Include' radio selected.")

    # 6. Expand "Other events"
    print("🔽 Expanding 'Other events'...")
    expand_btn = wait.until(EC.element_to_be_clickable((By.XPATH, "//span[contains(text(),'Other events')]/ancestor::div[contains(@class, 'mat-expansion-panel-header')]")))
    driver.execute_script("arguments[0].click();", expand_btn)
    time.sleep(1)

    # 7. Click "EM" checkbox
    print("☑️ Clicking EM checkbox...")
    em_checkbox = wait.until(EC.element_to_be_clickable((By.XPATH, "//span[contains(text(),'EM')]/preceding-sibling::div/input")))
    driver.execute_script("arguments[0].click();", em_checkbox)

    # 8. Click GET OASIS EVENTS
    print("📤 Clicking 'GET OASIS EVENTS'...")
    get_button = wait.until(EC.element_to_be_clickable((By.XPATH, "//button[.//span[contains(text(),'Get OASIS Events')]]")))
    driver.execute_script("arguments[0].click();", get_button)

    # 9. Wait for results
    wait.until(EC.presence_of_element_located((By.CLASS_NAME, "MuiTableBody-root")))
    print("✅ Results loaded.")

    # 10. Print rows
    print("\n📋 Retrieved Log Messages:\n")
    rows = driver.find_elements(By.CSS_SELECTOR, ".MuiTableBody-root > .MuiTableRow-root")
    for row in rows:
        cols = row.find_elements(By.TAG_NAME, "td")
        if len(cols) >= 3:
            print(f"{cols[0].text} | {cols[1].text} | {cols[2].text}")

except Exception as e:
    print(f"\n🚨 Error occurred: {e}")
    driver.save_screenshot("fatal_error.png")
    print("🛠️ Screenshot saved as 'fatal_error.png'")
finally:
    input("\n✅ Done. Press Enter to close browser...")
    driver.quit()

