import Icon from './svg/ovary.svg';

function OvaryIcon({ size }: { size: number }) {
  return (
    <div style={{ width: size }}>
      <img src={Icon} alt="ovary icon" />
    </div>
  );
}

export default OvaryIcon;
