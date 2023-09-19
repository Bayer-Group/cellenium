CREATE OR REPLACE FUNCTION current_user_groups()
    RETURNS text[]
    LANGUAGE plpython3u
    IMMUTABLE
AS
$$


"""
Getting a test token with Bambus:

conda activate bambus_user
from Bambus import azure_ad_auth
azure_ad_auth.current_accesstoken_azure_graphapi()
"""



from typing import List
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import io
import base64
from jose import jwt
import requests
import time

def get_session_token():
    import plpy
    r = plpy.execute("SELECT current_setting('postgraphile.auth_header_value', TRUE)::VARCHAR token")
    return [row for row in r][0]['token']


internal_employee_groups = [
    "BAG.U.2000.BAG.ALL_EMPLOYEES",
    "US_Bayer_All_Employees",
    "fi.u.us_bayer_all_employees"
]

# "Affiliates and joint laboratories" can have access to ABio datasets (Derek's email Oct 6th 2021)
abio_access_permitted_groups = [
    "BAG.U.2000.BAG.ALL_EMPLOYEES",
    "US_Bayer_All_Employees",
    "fi.u.us_bayer_all_employees",
    "e987bdec-2b99-480e-8c11-a8e7ef604a31",  # Azure AD partner-cellenium-askbio (Affiliate)
    "132b029b-bd15-4e2f-b8fb-5a0c321e4984"  # Azure AD partner-cellenium-broadinstitute (Joint Lab)
    # we don't include:
    # Recursion Pharma
]


def access_token_user_groups(access_token):
    global SD

    if not access_token:
        return ['UNAUTHENTICATED'];

    time_key = access_token + '_time'
    if time_key in SD:
        if time.monotonic() - SD[time_key] < 55 * 60:
            return SD[access_token]


    # TODO expire old entries...


    def _return_cached(access_token, groups):
        SD[access_token] = groups
        SD[time_key] = time.monotonic()
        return groups

    # plpy.warning(f'token: {access_token}')
    try:
        access_token_claims = jwt.get_unverified_claims(access_token)
        user_email = access_token_claims['unique_name']
        plpy.notice(f'user_email: {user_email}')
        ip_address = access_token_claims['ipaddr']

        headers = {
            'Authorization': f'Bearer {access_token}',
            'Host': 'graph.microsoft.com'
        }
        resp = requests.get('https://graph.microsoft.com/v1.0/me/memberOf', headers=headers)
        if resp.status_code >= 400:
            plpy.warning('response: %s' % resp)
            # log.error(resp.text)
            raise Exception('cannot get /me/memberOf')

        # getting the groups, with pagination
        groups = [d.get('mailNickname', d['displayName']) for d in resp.json()['value']]
        # also making group IDs available so that the membership in a partner group (which does not show its name in this API response) can be tied to study permissions
        groups.extend([d.get('id') for d in resp.json()['value']])
        while resp.json().get('@odata.nextLink') is not None:
            next_link = resp.json().get('@odata.nextLink')
            resp = requests.get(next_link, headers=headers)
            # 'id': d['id'],                          'b30e6aaf-b876-4a49-a376-5187bf958204'
            # 'displayName': d['displayName'],        'rdphsmi_sci'
            # 'mailNickname': d.get('mailNickname'),  'ph.a.rdphsmi_sci'
            groups.extend([d.get('mailNickname', d['displayName']) for d in resp.json()['value']])
            groups.extend([d.get('id') for d in resp.json()['value']])
        groups = [g for g in groups if g]
        groups.sort()

        # find department (only possible for bayer accounts)
        if not user_email.endswith('@bayer.com'):
            department_pseudogroups = ['partner group user']
        else:
            resp = requests.get('https://graph.microsoft.com/beta/me/profile/positions', headers=headers)
            if resp.status_code >= 400:
                # log.error('response: %s' % resp)
                # log.error(resp.text)
                raise Exception('cannot get /beta/me/profile/positions')
            department = None
            department_pseudogroups = []
            for e in resp.json()['value']:
                if e.get('detail') and e['detail'].get('company') and e['detail']['company'].get('department'):
                    department = e['detail']['company']['department']
            if department:
                # e.g. BAG-PH-RD-OIDT-DT-DTLB
                elements = department.split('-')
                department_pseudogroups = ['-'.join(elements[0:i + 1]) for i in range(0, (len(elements)))]

        final_groups = []
        final_groups.append(user_email)
        final_groups.extend(department_pseudogroups)
        final_groups.extend(groups)
        if len(set(groups) & set(internal_employee_groups)) > 0:
            final_groups.append('INTERNAL_EMPLOYEE')
        if len(set(groups) & set(abio_access_permitted_groups)) > 0:
            final_groups.append('ANALYTICAL_BIOSCIENCES_PERMITTED')
        plpy.notice('for %s (IP %s), found groups in Graph API: %s %s', user_email, ip_address, final_groups,
                    '(no-internal-employee-detected)' if 'INTERNAL_EMPLOYEE' not in final_groups else '')
        # log.info('for %s (IP %s), found groups in Graph API', user_email, ip_address)
        return _return_cached(access_token, final_groups)
    except Exception as e:
        plpy.warning(f'graph api {e}')
        return _return_cached(access_token, ['UNAUTHENTICATED'])


return access_token_user_groups(get_session_token())
$$;



CREATE OR REPLACE FUNCTION current_user_email()
    RETURNS text
    LANGUAGE plpython3u
    IMMUTABLE
AS
$$
from jose import jwt

def get_session_token():
    import plpy
    r = plpy.execute("SELECT current_setting('postgraphile.auth_header_value', TRUE)::VARCHAR token")
    return [row for row in r][0]['token']


def access_token_user_email(access_token):
    try:
        access_token_claims = jwt.get_unverified_claims(access_token)
        user_email = access_token_claims['unique_name']
        return user_email
    except Exception as e:
        plpy.warning(f'cannot retrieve email from JWT {e}')
    return 'anonymous'


return access_token_user_email(get_session_token())
$$;
